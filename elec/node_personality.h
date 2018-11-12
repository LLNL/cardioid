typedef struct _node_personality 
{
  Personality_t personality;
 
  /* global node id */
  uint32_t node_id;

  /* number of nodes */
  uint32_t nodes;

  /* node coordinates */
  uint32_t a, b, c, d, e;
  uint32_t coord[CR_NUM_DIMS];

  /* number of nodes in each dimensions */
  uint32_t a_nodes, b_nodes, c_nodes, d_nodes, e_nodes;
  uint32_t dim[CR_NUM_DIMS];

} node_personality;

