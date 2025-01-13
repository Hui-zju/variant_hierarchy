// Adjust SVG dimensions based on device size
const svg = d3.select("svg")
    .attr("width", window.innerWidth)
    .attr("height", window.innerHeight * 0.8)
    .call(d3.zoom().scaleExtent([0.1, 10]).on("zoom", (event) => g.attr("transform", event.transform))); // 缩放和平移 无用

// Margin and translation for tree rendering
const margin = { top: 20, right: 20, bottom: 20, left: 40 };
const width = +svg.attr("width") - margin.left - margin.right;
const height = +svg.attr("height") - margin.top - margin.bottom;
const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

// Tooltip setup
const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("position", "absolute")
    .style("visibility", "hidden")
    .style("background-color", "#f9f9f9")
    .style("border", "1px solid #ddd")
    .style("padding", "10px")
    .style("border-radius", "5px")
    .style("box-shadow", "0 2px 6px rgba(0,0,0,0.2)")
    .style("font-size", "12px");

// 搜索框元素
const geneSearchBox = document.getElementById("gene-search-box");
const geneSearchResults = document.getElementById("gene-search-results");

const searchBox = document.getElementById("search-box");
const searchResults = document.getElementById("search-results");
const nodeMap = new Map();

// Load CSV files and process data
Promise.all([
    d3.csv("node.csv"),
    d3.csv("edge.csv")
]).then(([nodes, edges]) => {
    nodes.forEach(node => {
        nodeMap.set(node["entity:ID"], { id: node["entity:ID"], name: node["name"], mutation: node["variant_info"] });
    });
    
    const edgeMap = edges.reduce((map, edge) => {
        if (!map[edge[":START_ID"]]) map[edge[":START_ID"]] = [];
        map[edge[":START_ID"]].push(edge[":END_ID"]);
        return map;
    }, {});


    function buildHierarchy(source, depth = 10) {
        const children = edgeMap[source] || [];
        depth = depth - 1;
        return depth > 1 ? {
            id: source,
            name: nodeMap.get(source)?.name || source,
            children: children.map(buildHierarchy),
            _children: null
        } : {
            id: source,
            name: nodeMap.get(source)?.name || source,
            children: null,
            _children: children.map(buildHierarchy)
        };
    }
    
    const initialTree = buildHierarchy("e2131", 2);  // 
    drawTree(initialTree);

    geneSearchBox.addEventListener("input", (event) => {
        const query = event.target.value.toLowerCase();
        if (query === "") {
            geneSearchResults.style.display = "none";
            return;
        }

        // const matches = nodes.filter(node => node.name.toLowerCase().includes(query));
        const matches = nodes
            .map(node => node.name.split(/[^a-zA-Z0-9]/)[0]) // Extract gene part
            .filter((gene, index, self) => gene.toLowerCase().includes(query) && self.indexOf(gene) === index);

        geneSearchResults.style.display = matches.length ? "block" : "none";
        geneSearchResults.innerHTML = matches.map(gene => `
            <div data-id="${gene}">${gene}</div>
        `).join("");

        Array.from(geneSearchResults.children).forEach(item => {
            item.addEventListener("click", () => {
                geneSearchBox.value = item.textContent;
                geneSearchBox.dataset.gene = item.textContent; // 保存选中的节点ID
                geneSearchResults.style.display = "none";

                const selectedGene = geneSearchBox.dataset.gene;
                const selectedNode = findNodeByGene(selectedGene)
                tree = buildHierarchy(selectedNode)
                drawTree(tree);
            });
        });
    });





    searchBox.addEventListener("input", (event) => {
        const query = event.target.value.toLowerCase();
        if (query === "") {
            searchResults.style.display = "none";
            return;
        }

        const matches = nodes.filter(node => node.name.toLowerCase().includes(query));
        searchResults.style.display = matches.length ? "block" : "none";
        searchResults.innerHTML = matches.map(node => `
            <div data-id="${node["entity:ID"]}">${node.name}</div>
        `).join("");

        Array.from(searchResults.children).forEach(item => {
            item.addEventListener("click", () => {
                searchBox.value = item.textContent;
                searchBox.dataset.selectedNode = item.dataset.id; // 保存选中的节点ID
                searchResults.style.display = "none";

                const selectedNodeId = searchBox.dataset.selectedNode;
                const selectedGene = searchBox.value.split(/[^a-zA-Z0-9]/)[0];
                const rootNode = findNodeByGene(selectedGene)
                tree = buildHierarchy(rootNode)
                drawTree(tree, selectedNodeId)
            });
        });
    });



    function findNodeByGene(gene) {
        const targetNodeName = `${gene} mutation`;
        const matchedNode = nodes.find(node => node["name"] === targetNodeName);
        if (matchedNode) {
            return matchedNode["entity:ID"];  // 返回找到的节点的 ID
        } else {
            return 'e0';  // 如果没有找到相应的节点，则返回根节点
        }
    }
});

function drawTree(data, focusedNodeId = data.id) {
    const tree = d3.tree().nodeSize([30, 200]);
    const root = d3.hierarchy(data, d => d.children || d._children);

    tree(root);
    const des = root.descendants();
    const findNode = (root, target) => root.descendants().find(d => d.data.id === target);
    let focus = findNode(root, focusedNodeId) || root;

    // 动态调整 SVG 大小
    const svgWidth = root.height * 200 + margin.left + margin.right;
    const svgHeight = root.descendants().length * 30 + margin.top + margin.bottom;
    // svg.attr("width", Math.max(svgWidth, window.innerWidth))
    //     .attr("height", Math.max(svgHeight, window.innerHeight * 0.8));

    const update = (focusedNode = root) => {
        const nodes = root.descendants();
        const links = root.links();

        // 清空画布并重新绘制
        g.selectAll("*").remove();

        // 绘制链接
        const link = g.selectAll(".link")
            .data(links)
            .join("path")
            .attr("class", "link")
            .attr("d", d3.linkHorizontal()
                .x(d => d.y)
                .y(d => d.x))
            .style("fill", "none")
            .style("stroke", "#ccc")
            .style("stroke-width", 1.5);

        // 绘制节点
        const node = g.selectAll(".node")
            .data(nodes)
            .join("g")
            .attr("class", "node")
            .attr("transform", d => `translate(${d.y},${d.x})`)
            .on("click", (event, d) => {
                if (d.children) {
                    d._children = d.children;
                    d.children = null; // 收起节点
                } else {
                    d.children = d._children;
                    d._children = null; // 展开节点
                }
                update(d);
            })
            .on("mouseover", (event, d) => {
                // Show the tooltip with mutation details
                const nodeData = nodeMap.get(d.data.id);
                
                // Parse the mutation string into an object
                let mutationData;
                try {
                    mutationData = JSON.parse(nodeData.mutation);
                } catch (error) {
                    mutationData = { displayName: "Unknown Mutation", type: "Unknown", source: "Unknown" };
                }

                tooltip.style("visibility", "visible")
                    .html(`
                        <strong>Mutation:</strong> ${mutationData.displayName} <br>
                        <strong>Gene:</strong> ${mutationData.reference1} <br>
                        ${mutationData.exon ? `<strong>Exon:</strong> ${mutationData.exon} <br>` : ""}
                        <strong>Type:</strong> ${mutationData.type} <br>
                        <strong>Source:</strong> ${mutationData.source.join(", ")}
                        
                    `)
                    .style("top", `${event.pageY + 10}px`)
                    .style("left", `${event.pageX + 10}px`);
            })
            .on("mouseout", () => {
                tooltip.style("visibility", "hidden");
            });

        node.append("circle")
            .attr("r", 5)
            .style("fill", d => (d.children || d._children) ? "#69b3a2" : "#ccc")
            .style("stroke", "steelblue")
            .style("stroke-width", 2);

        node.append("text")
            .attr("dy", 3)
            .attr("x", d => d.children ? -10 : 10)
            .style("text-anchor", d => d.children ? "end" : "start")
            .text(d => d.data.name)
            .style("font", "12px sans-serif");

        node.filter(d => d === focus)
        .select("circle")
        .style("fill", "red")
        .style("stroke", "darkred")
        .style("stroke-width", 3);

        node.filter(d => d === focus)
        .select("text")
        .style("font-weight", "bold")
        .style("font-size", "16px")
        .style("fill", "red");

        centerTree(focusedNode);
    };

    const centerTree = (focusedNode) => {
        const xOffset = height / 2 - focusedNode.x; // 竖直方向根据节点调整
        const yOffset = width / 3 - focusedNode.y; // 水平方向保持中心
        g.transition().duration(500).attr("transform", `translate(${yOffset},${xOffset})`);
    };

    update(focus);
}
