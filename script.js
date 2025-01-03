const nodeFile = 'node.csv';
const edgeFile = 'edge.csv';

// Load CSV data
Promise.all([d3.csv(nodeFile), d3.csv(edgeFile)]).then(([nodes, edges]) => {
  const nodeMap = new Map(nodes.map(node => [node["entity:ID"], node]));
  
  const treeData = buildTree(nodes, edges, nodeMap);
  drawTree(treeData);
});

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
    return {
      id,
      name: node.name,
      label: node[":LABEL"],
      children: (childrenMap[id] || []).map(createNode),
    };
  }

  return createNode(rootId);
}

function drawTree(data) {
  const width = 960;
  const height = 600;

  const treeLayout = d3.tree().size([height, width - 160]);
  const root = d3.hierarchy(data);

  const svg = d3.select('svg')
    .attr('viewBox', [0, 0, width, height]);

  const g = svg.append('g').attr('transform', 'translate(40,0)');

  treeLayout(root);

  const link = g.selectAll('.link')
    .data(root.links())
    .enter().append('path')
    .attr('class', 'link')
    .attr('d', d3.linkHorizontal()
      .x(d => d.y)
      .y(d => d.x));

  const node = g.selectAll('.node')
    .data(root.descendants())
    .enter().append('g')
    .attr('class', 'node')
    .attr('transform', d => `translate(${d.y},${d.x})`);

  node.append('circle')
    .attr('r', 5)
    .on('click', (event, d) => {
      d.children = d.children ? null : d._children;
      drawTree(data);
    });

  node.append('text')
    .attr('dy', 3)
    .attr('x', d => d.children ? -8 : 8)
    .style('text-anchor', d => d.children ? 'end' : 'start')
    .text(d => d.data.name);

  // Search functionality
  d3.select('#searchBox').on('input', function () {
    const searchTerm = this.value.toLowerCase();
    node.selectAll('circle')
      .classed('highlight', d => d.data.name.toLowerCase().includes(searchTerm));
  });
}