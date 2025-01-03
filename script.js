const nodeFile = 'node.csv';
const edgeFile = 'edge.csv';

// Load CSV data
Promise.all([d3.csv(nodeFile), d3.csv(edgeFile)]).then(([nodes, edges]) => {
  const data = {
    nodes: nodes.map(node => ({
      id: node["entity:ID"],
      name: node.name,
      label: node[":LABEL"]
    })),
    links: edges.map(edge => ({
      source: edge.source,
      target: edge.target,
      type: edge.contain
    }))
  };

  drawGraph(data);
});

function drawGraph(data) {
  const width = 800;
  const height = 600;

  const svg = d3.select('svg');
  const simulation = d3.forceSimulation(data.nodes)
    .force('link', d3.forceLink(data.links).id(d => d.id).distance(150))
    .force('charge', d3.forceManyBody().strength(-400))
    .force('center', d3.forceCenter(width / 2, height / 2));

  const link = svg.append('g')
    .selectAll('.link')
    .data(data.links)
    .enter().append('line')
    .attr('class', 'link');

  const node = svg.append('g')
    .selectAll('.node')
    .data(data.nodes)
    .enter().append('g')
    .attr('class', 'node')
    .call(d3.drag()
      .on('start', dragstarted)
      .on('drag', dragged)
      .on('end', dragended));

  node.append('circle')
    .attr('r', 10);

  node.append('text')
    .attr('dx', 12)
    .attr('dy', '.35em')
    .text(d => d.name);

  simulation.on('tick', () => {
    link
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);

    node.attr('transform', d => `translate(${d.x},${d.y})`);
  });

  // Search functionality
  d3.select('#searchBox').on('input', function () {
    const searchTerm = this.value.toLowerCase();
    node.selectAll('circle')
      .classed('highlight', d => d.name.toLowerCase().includes(searchTerm));
  });

  function dragstarted(event, d) {
    if (!event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
  }

  function dragged(event, d) {
    d.fx = event.x;
    d.fy = event.y;
  }

  function dragended(event, d) {
    if (!event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
  }
}
