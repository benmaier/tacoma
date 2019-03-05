class temporalNetworkFigure {

    constructor(
                config_file_name,
                container_id,
                scale = 1,
               )
    {
        var self = this;

        self.canvas = null;
        self.context = null;

        //var waiting_time = Math.pow(10,3-document.getElementById('log_dt').value);
        self.waiting_time = 20;
        self.timer = null;
        self.frame_skip = 1;
        self.paused = true;
        self.scale = scale;

        /*
            document.getElementById('log_dt').addEventListener('change', function() {
                console.log(+this.value);
                waiting_time = Math.pow(10,3-this.value);
                timer.stop();
                play();
            });
            document.getElementById('log_dt').addEventListener('input', function() {
                console.log(+this.value);
                waiting_time = Math.pow(10,3-this.value);
                console.log(Math.pow(10,-this.value));
                console.log(waiting_time);

                timer.stop();
                play();
            });


        // sociopatterns
        if (is_dtu){
            // DTU
            var link_distance = 20;
            var node_charge = -5;
            var file1 = "binned_dtu_1_week.taco";
            var file2 = "fw_binned_dtu_1_week.taco";
            var title_content_1 = "DTU Sep '14";
        } else {
            var link_distance = 20;
            var node_charge = -10;
            var file1 = "binned_ht09.taco";
            var file2 = "fw_binned_ht09.taco";
            var title_content_1 = "SocioPatterns HT09";
        }
        */
        self.min_waiting_time = 20;
        self.max_waiting_time = 1000;
        self.behavioral_change_scale = 0.6;
        self.waiting_time_scale = d3.scaleLog()
                                    .domain([0.001,self.behavioral_change_scale])
                                    .range([self.max_waiting_time, self.min_waiting_time]);
        self.frame_skip_scale = d3.scaleLinear()
                                    .domain([self.behavioral_change_scale,1])
                                    .range([1,20]);

        self.network_views = [];
        self.edge_views = [];
        self.edge_offsets_X = [];
        self.edge_offsets_Y = [];
        self.time = null;
        self.dt = null;
        self.plot_width;
        self.network_plot_height;

        self.mouseX = -1;
        self.mouseY = -1;
        self.global_it = 0;


        // load configuration data
        d3.queue()
            .defer(d3.json,config_file_name)
            .await(function(error, config) {

                var next_queue = d3.queue();

                var temporal_network_files = config.temporal_network_files;
                var edges_coordinate_files = config.edges_coordinate_files;

                if (config.hasOwnProperty('scale'))
                    self.scale = config.scale;

                self.plot_width = self.scale * config.plot_width;
                self.network_plot_height = self.scale * config.network_plot_height;

                self.global_width = self.plot_width * config.temporal_network_files.length;
                self.global_height = self.scale * (config.network_plot_height + config.edges_plot_height);

                for(var i = 0; i < temporal_network_files.length; i++)
                {
                    next_queue.defer(d3.json, "./"+temporal_network_files[i]);
                    next_queue.defer(d3.json, "./"+edges_coordinate_files[i]);
                }

                // load temporal_network data
                next_queue.await(function (error) {

                    self.canvas = d3.select('#'+container_id)
                      .append('canvas')
                      .attr("id","mainCanvas")
                      .attr('width', self.global_width)
                      .attr('height', self.global_height)
                    ;

                    self.context = self.canvas.node().getContext('2d');
                    self.context.fillStyle = "#fff"
                    self.context.fillRect(0,0,self.global_width,self.global_height);

                    self.global_it = config.start_it;

                    for (var i = 1; i < arguments.length; i++)
                    {
                        var network_count = Math.floor((i-1)/2);

                        var offset_X = self.plot_width * network_count;
                        if (i % 2 == 1)
                        {
                            var offset_Y = 0;
                            var tn_data = arguments[i];

                            if (i == 1)
                            {
                                self.time = tn_data.t;
                                self.dt = self.time[1] - self.time[0];
                            }

                            var this_tn = new temporalNetworkView(
                                            self.context,
                                            tn_data,
                                            config.titles[network_count],
                                            self.plot_width,
                                            self.network_plot_height,
                                            offset_X,
                                            offset_Y,
                                            self.scale * config.padding,
                                            self.scale * config.node_radius,
                                            self.scale * config.link_width,
                                            self.scale * config.node_edge_width,
                                            self.scale * config.font_size_in_px,
                                            self.scale * config.link_distance,
                                            self.scale * config.node_charge,
                                            1/self.scale * 0.1,
                                            config.start_it,
                                            config.d3_format_string
                                        );
                            self.network_views.push( this_tn );
                        }
                        else
                        {
                            var offset_Y = self.scale * config.network_plot_height;
                            var edge_data = arguments[i];
                            var this_ev = new temporalEdgesView(
                                                self.context,
                                                edge_data,
                                                self.plot_width,
                                                self.scale * config.edges_plot_height,
                                                offset_X,
                                                offset_Y,
                                                self.time[config.start_it],
                                                self.dt,
                                                self.scale * config.font_size_in_px / 1.5 + 1, // font size
                                                self.scale * config.padding,
                                                self.scale,
                                                self.scale * config.edge_line_width,
                                                config.d3_format_string
                                            );
                            self.edge_views.push( this_ev );

                            self.edge_offsets_X.push(offset_X);
                            self.edge_offsets_Y.push(offset_Y);
                        }
                    }
                    self.init_controls();
                    self.init_dragging();
                    self.play();
                });
            });
        }

        init_controls() {

            var self = this;

            self.canvas
                .on('mouseout', function () {
                    self.mouseX = -1;
                    self.mouseY = -1;
                    self.update_all_edge_draws();
                })
                .on('mousemove', function () {
                    self.mouseX = d3.event.layerX || d3.event.offsetX; 
                    self.mouseY = d3.event.layerY || d3.event.offsetY;
                    self.update_all_edge_draws();
                })
                .on('click', function () {
                    self.mouseX = d3.event.layerX || d3.event.offsetX; 
                    self.mouseY = d3.event.layerY || d3.event.offsetY;
                    var plot_id = Math.floor( self.mouseX / self.plot_width );

                    var this_ev = self.edge_views[plot_id];

                    if (self.mouseX > this_ev.padding + this_ev.offset_X && 
                        self.mouseX < this_ev.plot_width - 2*this_ev.padding + this_ev.offset_X &&
                        self.mouseY > this_ev.padding + this_ev.offset_Y && 
                        self.mouseY < this_ev.plot_height - 2*this_ev.padding + this_ev.offset_Y)
                    {
                        // get the current T marker position
                        var mouseX_in_original_coordinates = this_ev.inverse_fisheye_X(self.mouseX);
                        var this_mark_T = this_ev.xScale.invert(mouseX_in_original_coordinates);

                        // figure out whether or not the whole thing is paused and if it isn't, resume 
                        // playing afterwards
                        var is_currently_paused = self.paused;
                        if (!is_currently_paused)
                            self.pause();

                        //update the global time marker
                        self.global_it = Math.floor( this_mark_T / self.dt );
                        self.update_all_visualizations();

                        // resume playing if it was playing before
                        if (!is_currently_paused)
                            self.play();
                    }
                });
            ;/*
                */
        }

        pause()
        {
            if (!this.paused)
            {
                this.timer.stop();
                this.paused = true;
            }
        }

        play()
        {
            if (this.paused)
            {
                var self = this;
                self.paused = false;

                self.timer = d3.interval(function(){
                        self.update_frame();
                    }
                    ,self.waiting_time);
            }
        }

        update_frame()
        {
            var self = this;
            self.update_global_it();
            self.update_all_visualizations();
        }

        update_global_it(delta_it = null)
        {
            var self = this;

            if (delta_it === null) 
                delta_it = self.frame_skip;

            self.global_it += delta_it;

            if (self.global_it >= self.time.length)
                self.global_it = self.global_it % self.time.length;
            else if (self.global_it < 0)
                self.global_it += self.time.length;
        }

        update_aggregation_time(aggregation_time)
        {
            var self = this;

            self.network_views.forEach( function (tn) {
                tn.update_aggregation_time(aggregation_time);
            });
            self.edge_views.forEach( function (ev) {
                ev.update_aggregation_time(aggregation_time);
                /*
                ev.draw(  
                         self.mouseX, 
                         self.mouseY 
                       );
                */
            });

            if (self.paused)
            {
                self.update_all_visualizations();
            }
        }

        update_aggregation_time_scale(scale)
        {
            let self = this;
            let this_scale = d3.scalePow().exponent(3).domain([0.01,1]).range([self.dt, self.time[self.time.length-1]]);
            let aggregation_time = this_scale(scale);
            self.update_aggregation_time(aggregation_time);

            return aggregation_time;
        }

        update_animation_speed(scale)
        {
            var self = this;
            var abs_scale = Math.abs(scale);

            //console.log("scale =", scale);
            //console.log("abs_scale =", abs_scale);

            var is_currently_playing = ! self.paused;

            if (is_currently_playing)
                self.pause();

            // scale has to be a number between -1 and +1 
            if (abs_scale == 0)
            {
                if (!self.paused){
                    self.pause()
                }
                self.frame_skip = 0;
            }
            else 
            {
                var delta_it;
                if (abs_scale < self.behavioral_change_scale)
                {
                    self.waiting_time = self.waiting_time_scale(abs_scale);
                    delta_it = 1;
                }
                else
                {
                    self.waiting_time = self.min_waiting_time;
                    delta_it = Math.floor(self.frame_skip_scale(abs_scale));
                }
                var factor = 1;
                if (scale < 0)
                { 
                    factor = -1;

                    //console.log("scale < 0 ->", scale < 0);
                    //console.log("factor =", factor);
                }
                self.frame_skip = factor * delta_it;
            }

            //console.log("self.waiting_time =",self.waiting_time);
            //console.log("self.frame_skip =",self.frame_skip);

            if (is_currently_playing)
                self.play();

        }

        update_all_visualizations()
        {
            this.update_all_links();
            this.update_all_network_draws();
            this.update_all_edge_draws();
        }


        update_all_links() {
            var self = this;
            self.network_views.forEach( function (tn) {
                tn.update_links(self.global_it);
            });
        }


        update_all_network_draws() {
            var self = this;
            self.network_views.forEach( function (tn) {
                tn.update_draw();
            });
        }

        update_all_edge_draws() {
            var self = this;
            self.edge_views.forEach( function (ev) {
                ev.draw( self.time[self.global_it], 
                         self.mouseX, 
                         self.mouseY 
                       );
            });
        }

        init_dragging(){
            var canvas = document.querySelector("#mainCanvas");
            var self = this;
            console.log("init dragging");
            d3.select(canvas)
                .call(d3.drag()
                    .container(canvas)
                    .subject( function (){

                        var mouseX = d3.event.x; 
                        var mouseY = d3.event.y;
                        var plot_id = Math.floor( mouseX / self.plot_width );

                        self.current_drag_nv = self.network_views[plot_id];

                        if ( 
                            mouseY >= 0 && 
                            mouseY <= self.network_plot_height)
                        {
                            self.current_drag_nv = self.network_views[plot_id];
                            return self.current_drag_nv.simulation.find(d3.event.x - self.current_drag_nv.x0,
                                                                        d3.event.y - self.current_drag_nv.y0);
                        }
                        else
                        {
                            self.current_drag_nv = null;
                            return null;
                        }
                    })
                    .on("start", function () {
                        var simulation = self.current_drag_nv.simulation;
                        var n_view = self.current_drag_nv;
                        if (!d3.event.active) simulation.alphaTarget(n_view.alpha).restart();
                        d3.event.subject.fx = d3.event.subject.x;
                        d3.event.subject.fy = d3.event.subject.y;
                    })
                    .on("drag", function () {
                        var simulation = self.current_drag_nv.simulation;
                        var n_view = self.current_drag_nv;
                        var x = d3.event.x;
                        var y = d3.event.y;
                        if (x + n_view.x0 > n_view.offset_X + n_view.plot_width)
                            x = n_view.offset_X + n_view.plot_width - n_view.x0;
                        else if (x + n_view.x0 < n_view.offset_X)
                            x = n_view.offset_X - n_view.x0;
                        if (y + n_view.y0 > n_view.offset_Y + n_view.plot_height - 2*n_view.node_radius)
                            y = n_view.offset_Y + n_view.plot_height - n_view.y0 - 2*n_view.node_radius;
                        else if (y + n_view.y0 < n_view.offset_Y)
                            y = n_view.offset_Y - n_view.y0;

                        d3.event.subject.fx = x;
                        d3.event.subject.fy = y;
                    })
                    .on("end", function() {
                        var simulation = self.current_drag_nv.simulation;
                        var n_view = self.current_drag_nv;

                        if (!d3.event.active) simulation.alphaTarget(0);
                        d3.event.subject.fx = null;
                        d3.event.subject.fy = null;
                        self.current_drag_nv = null;
                    })
                );
        }

}
