class temporalNetworkView {

    constructor(
                context,
                network_data,
                title,
                plot_width, 
                plot_height, 
                offset_X,
                offset_Y,
                padding = 0,
                node_radius = 3,
                link_width = 2,
                node_edge_width = 1.3,
                font_size_in_px = 18,
                link_distance = 20,
                node_charge = -5,
                gravity_strength = 1,
                start_it = 0,
                d3_format_string = "r",
                aggregate_frames = 1,
               ) {

        this.context = context;
        this.data = network_data;
        this.title = title;
        this.plot_width = plot_width;
        this.plot_height = plot_height;
        this.offset_X = offset_X;
        this.offset_Y = offset_Y;
        this.node_radius = node_radius;
        this.link_width = link_width;
        this.node_edge_width = node_edge_width;
        this.font_size = font_size_in_px;
        this.link_distance = link_distance;
        this.node_charge = node_charge;
        this.padding = padding;
        this.d3_format_string = d3_format_string;

        this.show_title = true;

        this.len_t = this.data.t.length;
        this.dt = this.data.t[1] - this.data.t[0];
        this.it = start_it;
        this.alpha = 0.3;

        this.aggregate_frames = aggregate_frames;

        this.links = [];
        this.nodes = [];
        for(var n = 0; n < this.data.N; n++)
            this.nodes.push({"id": n});

        this.update_links(start_it);

        var self = this;
        this.x0 = this.plot_width / 2 + this.offset_X;
        this.y0 = this.plot_height / 2 + this.offset_Y;

        this.simulation = d3.forceSimulation(this.nodes)
            .force("charge", d3.forceManyBody().strength(this.node_charge))
            .force("link", d3.forceLink(this.links).distance(this.link_distance))
            .force("x", d3.forceX())
            .force("y", d3.forceY())
            .alphaTarget(this.alpha)
            .on("tick", function () { self.update_draw(); } );

        this.simulation.force("x").strength(gravity_strength);
        this.simulation.force("y").strength(gravity_strength);
       
    }

    update_aggregation_time(aggregation_time) {
        this.aggregate_frames = Math.floor( aggregation_time / this.dt );
    }

    update_links(it) {
        this.it = it;

        let these_links = this.links;
        let these_nodes = this.nodes;
        these_links.length = 0;

        let graph = {};
        these_nodes.forEach( function(n) {
            graph[n.id] = [];
        });


        for(var this_it = it; (this_it < it + this.aggregate_frames) && (this_it < this.data.t.length); this_it++)
        {
            let edges = this.data.edges[this_it];
            edges.forEach(function(edge){
                let source = edge[0];
                let target = edge[1];
                if (graph[source].indexOf(target) < 0)
                {
                    these_links.push({source: these_nodes[source], 
                                      target: these_nodes[target]});
                    graph[source].push(target);
                }
            });
        }

    }

    update_draw() {
        this.context.clearRect(this.offset_X,
                               this.offset_Y, 
                               this.plot_width, 
                               this.plot_height);
        this.context.fillStyle = "#fff"
        this.context.fillRect(this.offset_X,
                               this.offset_Y, 
                               this.plot_width, 
                               this.plot_height);
        this.context.setLineDash([]);
        
        this.context.font = String(this.font_size) + "px Arial";

        this.context.save()
        this.context.textAlign = "right";
        this.context.fillStyle = '#000';
        var text = d3.format(this.d3_format_string)(this.data.t[this.it]) + this.data.time_unit;
        this.context.fillText(text,
                                   this.plot_width + this.offset_X - 2*this.padding,
                                   this.plot_height + this.offset_Y - 2.2*this.font_size);
                                    /*
                                  this.offset_X + this.plot_width - 2* this.padding,
                                  this.offset_Y + this.font_size);
                                  */

        this.context.restore();
        if (this.show_title)
        {
            this.context.save()
            this.context.textAlign = "right";
            this.context.fillStyle = '#000';
            this.context.fillText(this.title,
                                   this.plot_width + this.offset_X - 2*this.padding,
                                   this.plot_height + this.offset_Y - this.font_size);
            this.context.restore();
        }

        this.context.fillStyle = '#000'
        this.context.strokeStyle = '#fff'
        this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
        this.context.lineWidth = this.link_width;
        this.context.beginPath();

        let this_context = this.context;
        let self = this;

        this.links.forEach( function (link) {
                var x1 = link.source.x + self.x0;
                var y1 = link.source.y + self.y0;
                var x2 = link.target.x + self.x0;
                var y2 = link.target.y + self.y0;
                this_context.moveTo(x1,y1);
                this_context.lineTo(x2,y2);
        });

        this.context.stroke();

        this.context.beginPath();
        this.context.lineWidth = this.node_edge_width;
        this.context.fillStyle = '#000000';
        this.context.strokeStyle = '#ffffff';
        this.nodes.forEach( function (node) {
            this_context.moveTo(node.x + self.node_radius + self.x0, node.y + self.y0);
            this_context.arc(node.x + self.x0 ,node.y + self.y0,self.node_radius,0,2*Math.PI);
        });
        this.context.closePath();
        this.context.stroke();
        this.context.fill();

        this.simulation.nodes(this.nodes);
        this.simulation.force("link").links(this.links);
        this.simulation.alpha(this.alpha).restart();
        
    }
}
