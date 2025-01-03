const nodeFile = 'node.csv';
const edgeFile = 'edge.csv';

// Load CSV data
Promise.all([d3.csv(nodeFile), d3.csv(edgeFile)]).then(([nodes, edges]) => {
  const nodeMap = new Map(nodes.map(node => [node["entity:ID"], node]));
  const treeData = buildTree(nodes, edges, nodeMap);
  drawTree(treeData);
});

// Build hierarchical tree data
function buildTree(nodes, edges, nodeMap) {
  const rootId = nodes[0]["entity:ID"];
  const childrenMap = edges.reduce((map, edge) => {
    const parent = edge.source;
    const child = edge.target;
    if (!map[parent]) map[parent] = [];
    map[parent].push(child);
    return map;
  }, {});

  function createNode(id) {
    const node = nodeMap.get(id);
    const children = (childrenMap[id] || []).map(createNode);
    return {
      id,
      name: node.name,
      label: node[":LABEL"],
      children: children.length > 0 ? null : children, // Initially collapsed
      _children: children.length > 0 ? children : null // Store collapsed children
    };
  }

  return createNode(rootId);
}

// Draw the tree
function drawTree(data) {
  const width = 960;
  const height = 600;

  const treeLayout = d3.tree().size([height, width - 160]);
  const root = d3.hierarchy(data, d => d.children || d._children);
  treeLayout(root);

  const svg = d3.select('svg')
    .attr('viewBox', [0, 0, width, height])
    .selectAll('*')
    .remove(); // Clear previous tree

  const g = d3.select('svg')
    .append('g')
    .attr('transform', 'translate(40,0)');

  // Links
  g.selectAll('.link')
    .data(root.links())
    .enter().append('path')
    .attr('class', 'link')
    .attr('d', d3.linkHorizontal()
      .x(d => d.y)
      .y(d => d.x))
    .style('fill', 'none')
    .style('stroke', '#ccc')
    .style('stroke-width', 1.5);

  // Nodes
  const node = g.selectAll('.node')
    .data(root.descendants())
    .enter().append('g')
    .attr('class', 'node')
    .attr('transform', d => `translate(${d.y},${d.x})`);

  // Node circles
  node.append('circle')
    .attr('r', 5)
    .style('fill', d => d._children ? '#69b3a2' : '#fff')
    .style('stroke', '#555')
    .style('stroke-width', 1.5)
    .on('click', (event, d) => {
      toggleNode(d);
      drawTree(data); // Re-render the tree
    });

  // Node labels
  node.append('text')
    .attr('dy', 3)
    .attr('x', d => d.children || d._children ? -8 : 8)
    .style('text-anchor', d => d.children || d._children ? 'end' : 'start')
    .text(d => d.data.name);
}

// Toggle node expand/collapse
function toggleNode(d) {
  if (d.children) {
    d._children = d.children; // Collapse children
    d.children = null;
  } else {
    d.children = d._children; // Expand children
    d._children = null;
  }
}
