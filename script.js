// Load CSV files and process data
Promise.all([
    d3.csv("node.csv"),
    d3.csv("edge.csv")
]).then(([nodes, edges]) => {
    const nodeMap = new Map();
    nodes.forEach(node => {
        nodeMap.set(node.entity, { id: node.entity, name: node.name });
    });

    const root = { id: "root", name: "Root", children: [] };
    const edgeMap = edges.reduce((map, edge) => {
        if (!map[edge.source]) map[edge.source] = [];
        map[edge.source].push(edge.target);
        return map;
    }, {});

    // Build hierarchy
    function buildHierarchy(source) {
        const children = edgeMap[source] || [];
        return {
            id: source,
            name: nodeMap.get(source)?.name || source,
            children: children.map(buildHierarchy)
        };
    }
    root.children = buildHierarchy("e0").children;

    // Visualize with D3.js
    drawTree(root);
});

function drawTree(data) {
    const width = 960;
    const height = 600;
    const svg = d3.select("svg");
    const g = svg.append("g").attr("transform", "translate(40,0)");

    const tree = d3.tree().size([height, width - 160]);
    const root = d3.hierarchy(data);

    tree(root);

    const link = g.selectAll(".link")
        .data(root.links())
        .enter().append("path")
        .attr("class", "link")
        .attr("d", d3.linkHorizontal()
            .x(d => d.y)
            .y(d => d.x));

    const node = g.selectAll(".node")
        .data(root.descendants())
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", d => `translate(${d.y},${d.x})`);

    node.append("circle")
        .attr("r", 5);

    node.append("text")
        .attr("dy", 3)
        .attr("x", d => d.children ? -10 : 10)
        .style("text-anchor", d => d.children ? "end" : "start")
        .text(d => d.data.name);
}