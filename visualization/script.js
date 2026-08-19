const MAX_VISIBLE_SIBLINGS = 5;
const NODE_CSV_PATH = "../data/graph/node_with_info.csv";
const EDGE_CSV_PATH = "../data/graph/edge.csv";
const PREFERRED_ROOT_ID = "e2131";

const margin = { top: 24, right: 48, bottom: 24, left: 48 };
const svg = d3.select("#tree-svg");
const g = svg.append("g");
const zoom = d3.zoom()
    .scaleExtent([0.1, 10])
    .on("zoom", (event) => g.attr("transform", event.transform));

svg.call(zoom);

const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("visibility", "hidden");

const geneSearchBox = document.getElementById("gene-search-box");
const geneSearchResults = document.getElementById("gene-search-results");
const geneSearchWrapper = document.getElementById("gene-search-wrapper");
const searchBox = document.getElementById("search-box");
const searchResults = document.getElementById("search-results");
const nodeSearchWrapper = document.getElementById("node-search-wrapper");
const viewModeSelect = document.getElementById("view-mode");
const pathControl = document.getElementById("path-control");
const pathSelect = document.getElementById("path-select");
const resetViewButton = document.getElementById("reset-view");
const statusBar = document.getElementById("status-bar");

const state = {
    nodes: [],
    nodeMap: new Map(),
    childrenMap: new Map(),
    parentMap: new Map(),
    genes: [],
    selectedGene: "",
    selectedNodeId: null,
    viewMode: "gene",
    pathIndex: 0,
    ancestorPaths: [],
    expandedSiblingGroups: new Set(),
    collapsedNodes: new Set(),
    focusNodeId: null,
    viewportAnchor: null
};

Promise.all([
    d3.csv(NODE_CSV_PATH),
    d3.csv(EDGE_CSV_PATH)
]).then(([nodes, edges]) => {
    state.nodes = nodes;

    nodes.forEach(node => {
        const id = node["entity:ID"];
        const variant = parseVariant(node.variant_info);
        state.nodeMap.set(id, {
            id,
            name: node.name || id,
            mutation: node.variant_info,
            variant,
            gene: getGeneFromNode(node, variant)
        });
    });

    edges.forEach(edge => {
        addToMapList(state.childrenMap, edge[":START_ID"], edge[":END_ID"]);
        addToMapList(state.parentMap, edge[":END_ID"], edge[":START_ID"]);
    });

    state.genes = getGenes();
    bindControls();

    const initialRootId = pickInitialRootId();
    const initialNode = state.nodeMap.get(initialRootId);
    state.selectedGene = initialNode?.gene || getGeneFromName(initialNode?.name || "") || "";
    state.focusNodeId = initialRootId;
    geneSearchBox.value = state.selectedGene;
    renderCurrentView();
}).catch(error => {
    statusBar.textContent = `Failed to load data: ${error.message}`;
});

function bindControls() {
    geneSearchBox.addEventListener("input", () => {
        renderSearchResults(geneSearchResults, getGeneMatches(geneSearchBox.value), gene => {
            geneSearchBox.value = gene;
            searchBox.value = "";
            state.selectedGene = gene;
            state.selectedNodeId = null;
            state.viewMode = "gene";
            viewModeSelect.value = "gene";
            resetTransientViewState();
            state.focusNodeId = findNodeByGene(gene);
            renderCurrentView();
        });
    });

    searchBox.addEventListener("input", () => {
        const query = searchBox.value.trim().toLowerCase();
        const matches = query
            ? state.nodes.filter(node => (node.name || "").toLowerCase().includes(query)).slice(0, 80)
            : [];

        renderSearchResults(searchResults, matches, node => {
            const nodeId = node["entity:ID"];
            const nodeData = state.nodeMap.get(nodeId);

            searchBox.value = node.name;
            state.selectedNodeId = nodeId;
            state.selectedGene = nodeData?.gene || getGeneFromName(node.name);
            geneSearchBox.value = state.selectedGene;
            state.viewMode = "node";
            viewModeSelect.value = "node";
            state.focusNodeId = nodeId;
            resetTransientViewState();
            renderCurrentView();
        }, node => node.name);
    });

    viewModeSelect.addEventListener("change", () => {
        state.viewMode = viewModeSelect.value;
        resetTransientViewState();
        state.focusNodeId = state.selectedNodeId || findNodeByGene(state.selectedGene);
        geneSearchResults.style.display = "none";
        searchResults.style.display = "none";
        renderCurrentView();
    });

    pathSelect.addEventListener("change", () => {
        state.pathIndex = Number(pathSelect.value) || 0;
        state.expandedSiblingGroups.clear();
        state.focusNodeId = state.selectedNodeId;
        renderCurrentView();
    });

    resetViewButton.addEventListener("click", () => {
        state.expandedSiblingGroups.clear();
        state.collapsedNodes.clear();
        state.focusNodeId = state.selectedNodeId || findNodeByGene(state.selectedGene);
        renderCurrentView();
    });

    window.addEventListener("resize", renderCurrentView);
    document.addEventListener("click", event => {
        if (!event.target.closest(".search-wrapper")) {
            geneSearchResults.style.display = "none";
            searchResults.style.display = "none";
        }
    });
}

function renderCurrentView() {
    updateVisibleControls();
    setSvgSize();

    const selectedNode = state.selectedNodeId ? state.nodeMap.get(state.selectedNodeId) : null;
    const defaultRootId = pickInitialRootId();
    const geneRootId = state.selectedGene ? findNodeByGene(state.selectedGene) : defaultRootId;
    const effectiveRootId = geneRootId || defaultRootId;

    let treeData;
    let focusId;

    if (state.viewMode === "node" && !selectedNode) {
        state.ancestorPaths = [];
        state.pathIndex = 0;
        updatePathControl();
        renderEmptyView("Node hierarchy: search for a variant node to visualize.");
        return;
    }

    if (state.viewMode === "node") {
        state.ancestorPaths = collectAncestorPaths(state.selectedNodeId, effectiveRootId);
        state.pathIndex = clamp(state.pathIndex, 0, Math.max(state.ancestorPaths.length - 1, 0));
        treeData = buildNodeHierarchy(state.ancestorPaths[state.pathIndex] || [state.selectedNodeId]);
        focusId = state.focusNodeId || state.selectedNodeId;
    } else {
        const focusPath = selectedNode
            ? findDescendantPath(effectiveRootId, state.selectedNodeId)
            : [effectiveRootId];
        const pathChildByParent = makePathChildMap(focusPath);

        treeData = buildDescendantHierarchy(effectiveRootId, pathChildByParent);
        focusId = state.focusNodeId || (selectedNode && focusPath.length ? state.selectedNodeId : effectiveRootId);
        state.ancestorPaths = [];
        state.pathIndex = 0;
    }

    if (!treeContainsId(treeData, focusId)) {
        focusId = selectedNode && treeContainsId(treeData, state.selectedNodeId) ? state.selectedNodeId : treeData.id;
        state.focusNodeId = focusId;
    }

    updatePathControl();
    updateStatus(treeData, focusId);
    drawTree(treeData, focusId);
}

function updateVisibleControls() {
    const isNodeMode = state.viewMode === "node";

    geneSearchWrapper.hidden = isNodeMode;
    nodeSearchWrapper.hidden = !isNodeMode;
    pathControl.hidden = !isNodeMode;
}

function buildNodeHierarchy(path) {
    const selectedId = path[path.length - 1];
    const selectedParentId = path.length > 1 ? path[path.length - 2] : null;

    function buildPathNode(index) {
        const id = path[index];
        const node = createTreeNode(id);

        if (state.collapsedNodes.has(id)) {
            node.children = null;
            node.hiddenChildCount = getChildIds(id).length;
            return node;
        }

        if (id === selectedId) {
            const descendant = buildDescendantHierarchy(id, new Map());
            node.children = descendant.children;
            node.hiddenChildCount = descendant.hiddenChildCount;
            return node;
        }

        const nextId = path[index + 1];

        if (id !== selectedParentId) {
            node.children = [buildPathNode(index + 1)];
            return node;
        }

        const childIds = getChildIds(id);
        const visibleChildren = getVisibleChildIds(id, childIds, selectedId);
        node.children = visibleChildren.map(childId => {
            if (childId === nextId) {
                return buildPathNode(index + 1);
            }

            const sibling = createTreeNode(childId);
            sibling.contextOnly = true;
            return sibling;
        });

        appendMoreNode(node, id, childIds, visibleChildren);
        return node;
    }

    return buildPathNode(0);
}

function buildDescendantHierarchy(id, pathChildByParent, visited = new Set()) {
    const node = createTreeNode(id);
    const childIds = getChildIds(id);
    node.hasChildren = childIds.length > 0;

    if (visited.has(id) || state.collapsedNodes.has(id)) {
        node.children = null;
        node.hiddenChildCount = childIds.length;
        return node;
    }

    const nextVisited = new Set(visited);
    nextVisited.add(id);

    const forcedChildId = pathChildByParent.get(id);
    const visibleChildren = getVisibleChildIds(id, childIds, forcedChildId);
    node.children = visibleChildren.map(childId => buildDescendantHierarchy(childId, pathChildByParent, nextVisited));
    appendMoreNode(node, id, childIds, visibleChildren);

    return node;
}

function drawTree(data, focusedNodeId = data.id) {
    const width = Number(svg.attr("width"));
    const height = Number(svg.attr("height"));
    const tree = d3.tree().nodeSize([34, 220]);
    const root = d3.hierarchy(data, d => d.children);

    tree(root);

    const nodes = root.descendants();
    const links = root.links();
    const focus = nodes.find(d => d.data.id === focusedNodeId) || root;

    svg.attr("viewBox", `0 0 ${width} ${height}`);
    g.selectAll("*").remove();

    g.selectAll(".link")
        .data(links)
        .join("path")
        .attr("class", "link")
        .attr("d", d3.linkHorizontal()
            .x(d => d.y)
            .y(d => d.x));

    const node = g.selectAll(".node")
        .data(nodes)
        .join("g")
        .attr("class", d => `node${d.data.isMore ? " node-more" : ""}${d.data.contextOnly ? " node-context" : ""}`)
        .attr("transform", d => `translate(${d.y},${d.x})`)
        .on("click", (event, d) => {
            event.stopPropagation();

            if (d.data.isMore) {
                captureViewportAnchor(d.parent || d);
                state.expandedSiblingGroups.add(d.data.parentId);
                renderCurrentView();
                return;
            }

            if (!d.data.hasChildren && !d.data.hiddenChildCount) {
                return;
            }

            captureViewportAnchor(d);

            if (state.collapsedNodes.has(d.data.id)) {
                state.collapsedNodes.delete(d.data.id);
            } else {
                state.collapsedNodes.add(d.data.id);
            }

            state.focusNodeId = d.data.id;
            renderCurrentView();
        })
        .on("mouseover", showTooltip)
        .on("mousemove", moveTooltip)
        .on("mouseout", hideTooltip);

    node.append("circle")
        .attr("r", d => d.data.isMore ? 18 : 16)
        .attr("class", "hit-circle");

    node.append("circle")
        .attr("r", d => d.data.isMore ? 7 : 5)
        .attr("class", d => {
            if (d.data.isMore) return "more-circle";
            if (d.data.id === focusedNodeId) return "focus-circle";
            if (d.data.hasChildren || d.data.hiddenChildCount) return "branch-circle";
            return "leaf-circle";
        });

    node.append("text")
        .attr("dy", 4)
        .attr("x", 12)
        .attr("text-anchor", "start")
        .attr("class", d => d.data.id === focusedNodeId ? "focus-label" : "")
        .text(d => d.data.name);

    positionView(focus, nodes, width, height);
}

function renderEmptyView(message) {
    const width = Number(svg.attr("width"));
    const height = Number(svg.attr("height"));

    svg.attr("viewBox", `0 0 ${width} ${height}`);
    g.selectAll("*").remove();
    tooltip.style("visibility", "hidden");
    statusBar.textContent = message;

    svg.interrupt().call(zoom.transform, d3.zoomIdentity);
}

function positionView(focusedNode, nodes, width, height) {
    if (state.viewportAnchor) {
        const anchor = state.viewportAnchor;
        const anchorNode = nodes.find(node => node.data.id === anchor.id);
        state.viewportAnchor = null;

        if (anchorNode) {
            const transform = d3.zoomIdentity
                .translate(anchor.screenX - anchorNode.y * anchor.scale, anchor.screenY - anchorNode.x * anchor.scale)
                .scale(anchor.scale);

            svg.interrupt().call(zoom.transform, transform);
            return;
        }
    }

    const scale = 1;
    const xOffset = width / 3 - focusedNode.y * scale;
    const yOffset = height / 2 - focusedNode.x * scale;
    const transform = d3.zoomIdentity.translate(xOffset, yOffset).scale(scale);

    svg.transition().duration(300).call(zoom.transform, transform);
}

function captureViewportAnchor(d) {
    const transform = d3.zoomTransform(svg.node());

    state.viewportAnchor = {
        id: d.data.id,
        screenX: transform.applyX(d.y),
        screenY: transform.applyY(d.x),
        scale: transform.k
    };
}

function renderSearchResults(container, items, onPick, getLabel = item => item) {
    if (!items.length) {
        container.style.display = "none";
        container.innerHTML = "";
        return;
    }

    container.innerHTML = items.map((item, index) => {
        const label = escapeHtml(getLabel(item));
        return `<button type="button" data-index="${index}">${label}</button>`;
    }).join("");
    container.style.display = "block";

    Array.from(container.children).forEach(button => {
        button.addEventListener("click", () => {
            container.style.display = "none";
            onPick(items[Number(button.dataset.index)]);
        });
    });
}

function updatePathControl() {
    if (state.viewMode !== "node") {
        pathControl.hidden = true;
        pathSelect.disabled = true;
        pathSelect.innerHTML = `<option>Only available in node mode</option>`;
        return;
    }

    pathControl.hidden = false;

    if (!state.selectedNodeId) {
        pathSelect.disabled = true;
        pathSelect.innerHTML = `<option>Search a node first</option>`;
        return;
    }

    if (state.ancestorPaths.length <= 1) {
        const path = state.ancestorPaths[0] || [state.selectedNodeId];
        const label = path.map(id => state.nodeMap.get(id)?.name || id).join(" > ");
        pathSelect.disabled = true;
        pathSelect.innerHTML = `<option>${escapeHtml(label || "Single parent path")}</option>`;
        return;
    }

    pathSelect.disabled = false;
    pathSelect.innerHTML = state.ancestorPaths.map((path, index) => {
        const label = path.map(id => state.nodeMap.get(id)?.name || id).join(" > ");
        return `<option value="${index}">${escapeHtml(label)}</option>`;
    }).join("");
    pathSelect.value = String(state.pathIndex);
}

function updateStatus(treeData, focusId) {
    const focusName = state.nodeMap.get(focusId)?.name || focusId;
    const rootName = treeData?.name || "";
    const modeText = state.viewMode === "node" ? "Node hierarchy" : "Gene hierarchy";
    const pathText = state.ancestorPaths.length > 1
        ? ` Path ${state.pathIndex + 1}/${state.ancestorPaths.length}.`
        : "";

    statusBar.textContent = `${modeText}: ${rootName}. Focus: ${focusName}.${pathText}`;
}

function showTooltip(event, d) {
    if (d.data.isMore) {
        tooltip.style("visibility", "visible")
            .html(`<strong>${escapeHtml(d.data.name)}</strong><br>Click to show hidden siblings.`);
        moveTooltip(event);
        return;
    }

    const nodeData = state.nodeMap.get(d.data.id);
    const mutationData = nodeData?.variant || {};
    const source = Array.isArray(mutationData.source)
        ? mutationData.source.join(", ")
        : mutationData.source || "Unknown";
    const exon = getVariantExon(mutationData);

    tooltip.style("visibility", "visible")
        .html(`
            <strong>Mutation:</strong> ${escapeHtml(getVariantDisplayName(mutationData, nodeData?.name))}<br>
            <strong>Gene:</strong> ${escapeHtml(getGeneFromVariant(mutationData) || nodeData?.gene || "Unknown")}<br>
            ${exon ? `<strong>Exon:</strong> ${escapeHtml(exon)}<br>` : ""}
            <strong>Type:</strong> ${escapeHtml(getVariantType(mutationData))}<br>
            <strong>Source:</strong> ${escapeHtml(source)}
        `);
    moveTooltip(event);
}

function moveTooltip(event) {
    tooltip
        .style("top", `${event.pageY + 12}px`)
        .style("left", `${event.pageX + 12}px`);
}

function hideTooltip() {
    tooltip.style("visibility", "hidden");
}

function getVisibleChildIds(parentId, childIds, forcedChildId) {
    const ordered = prioritizeVisibleChildren(childIds, forcedChildId);

    if (childIds.length <= MAX_VISIBLE_SIBLINGS || state.expandedSiblingGroups.has(parentId)) {
        return ordered;
    }

    const visible = ordered.slice(0, MAX_VISIBLE_SIBLINGS);

    return [...new Set(visible)];
}

function prioritizeVisibleChildren(childIds, forcedChildId) {
    const childIndex = new Map(childIds.map((id, index) => [id, index]));

    return childIds.slice().sort((a, b) => {
        if (forcedChildId) {
            if (a === forcedChildId) return -1;
            if (b === forcedChildId) return 1;
        }

        const aHasChildren = getChildIds(a).length > 0;
        const bHasChildren = getChildIds(b).length > 0;

        if (aHasChildren !== bHasChildren) {
            return aHasChildren ? -1 : 1;
        }

        return childIndex.get(a) - childIndex.get(b);
    });
}

function appendMoreNode(node, parentId, childIds, visibleChildren) {
    if (childIds.length <= visibleChildren.length || state.expandedSiblingGroups.has(parentId)) {
        return;
    }

    const hiddenCount = childIds.length - visibleChildren.length;
    if (!node.children) {
        node.children = [];
    }

    node.children.push({
        id: `__more__${parentId}`,
        name: `Show ${hiddenCount} more siblings`,
        parentId,
        isMore: true,
        hasChildren: false,
        children: null
    });
}

function collectAncestorPaths(targetId, preferredRootId) {
    const paths = [];
    const stack = [{ id: targetId, path: [targetId], seen: new Set([targetId]) }];

    while (stack.length && paths.length < 80) {
        const current = stack.pop();
        const parents = getParentIds(current.id);

        if (current.id === preferredRootId || !parents.length) {
            paths.push(current.path.slice().reverse());
            continue;
        }

        parents.forEach(parentId => {
            if (!current.seen.has(parentId)) {
                stack.push({
                    id: parentId,
                    path: [...current.path, parentId],
                    seen: new Set([...current.seen, parentId])
                });
            }
        });
    }

    return paths.sort((a, b) => {
        const aStartsPreferred = a[0] === preferredRootId ? 0 : 1;
        const bStartsPreferred = b[0] === preferredRootId ? 0 : 1;
        return aStartsPreferred - bStartsPreferred || a.length - b.length;
    });
}

function findDescendantPath(rootId, targetId) {
    if (!rootId || !targetId) {
        return [];
    }

    const stack = [{ id: rootId, path: [rootId], seen: new Set([rootId]) }];

    while (stack.length) {
        const current = stack.pop();
        if (current.id === targetId) {
            return current.path;
        }

        getChildIds(current.id).forEach(childId => {
            if (!current.seen.has(childId)) {
                stack.push({
                    id: childId,
                    path: [...current.path, childId],
                    seen: new Set([...current.seen, childId])
                });
            }
        });
    }

    return [];
}

function makePathChildMap(path) {
    const pathChildByParent = new Map();

    for (let index = 0; index < path.length - 1; index += 1) {
        pathChildByParent.set(path[index], path[index + 1]);
    }

    return pathChildByParent;
}

function createTreeNode(id) {
    const nodeData = state.nodeMap.get(id);
    const childIds = getChildIds(id);

    return {
        id,
        name: nodeData?.name || id,
        hasChildren: childIds.length > 0,
        children: null
    };
}

function treeContainsId(node, id) {
    if (!node || !id) {
        return false;
    }

    if (node.id === id) {
        return true;
    }

    return (node.children || []).some(child => treeContainsId(child, id));
}

function findNodeByGene(gene) {
    const exactName = `${gene} mutation`.toLowerCase();
    const exactNode = state.nodes.find(node => (node.name || "").toLowerCase() === exactName);

    if (exactNode) {
        return exactNode["entity:ID"];
    }

    const geneNode = state.nodes.find(node => {
        const data = state.nodeMap.get(node["entity:ID"]);
        return data?.gene?.toLowerCase() === gene.toLowerCase() && getVariantType(data.variant) === "mutation";
    });

    return geneNode?.["entity:ID"] || pickInitialRootId();
}

function getGenes() {
    return [...new Set([...state.nodeMap.values()].map(node => node.gene).filter(Boolean))]
        .sort((a, b) => a.localeCompare(b));
}

function getGeneMatches(query) {
    const normalized = query.trim().toLowerCase();
    if (!normalized) {
        return [];
    }

    return state.genes
        .filter(gene => gene.toLowerCase().includes(normalized))
        .slice(0, 80);
}

function getGeneFromNode(node, variant) {
    return getGeneFromVariant(variant) || getGeneFromName(node.name || "");
}

function getGeneFromName(name) {
    return (name.match(/^[a-zA-Z0-9]+/) || [""])[0];
}

function getGeneFromVariant(variant) {
    return variant?.location?.gene || "";
}

function getVariantDisplayName(variant, fallbackName = "Unknown") {
    return variant?.display_name || fallbackName || "Unknown";
}

function getVariantType(variant) {
    return variant?.alteration?.type || variant?.variant_class || "Unknown";
}

function getVariantExon(variant) {
    return variant?.location?.exon || "";
}

function pickInitialRootId() {
    if (state.nodeMap.has(PREFERRED_ROOT_ID)) {
        return PREFERRED_ROOT_ID;
    }

    const root = state.nodes.find(node => !state.parentMap.has(node["entity:ID"]));
    return root?.["entity:ID"] || state.nodes[0]?.["entity:ID"] || "";
}

function parseVariant(value) {
    if (!value) {
        return {};
    }

    try {
        return JSON.parse(value);
    } catch (error) {
        return {};
    }
}

function getChildIds(id) {
    return state.childrenMap.get(id) || [];
}

function getParentIds(id) {
    return state.parentMap.get(id) || [];
}

function addToMapList(map, key, value) {
    if (!map.has(key)) {
        map.set(key, []);
    }

    map.get(key).push(value);
}

function resetTransientViewState() {
    state.pathIndex = 0;
    state.ancestorPaths = [];
    state.expandedSiblingGroups.clear();
    state.collapsedNodes.clear();
}

function setSvgSize() {
    const controlsHeight = document.getElementById("top-panel").offsetHeight;
    const width = Math.max(document.getElementById("visualization").clientWidth, 320);
    const height = Math.max(window.innerHeight - controlsHeight - 32, 520);

    svg.attr("width", width).attr("height", height);
}

function clamp(value, min, max) {
    return Math.min(Math.max(value, min), max);
}

function escapeHtml(value) {
    return String(value ?? "")
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#039;");
}
