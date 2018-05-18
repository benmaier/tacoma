class temporalNetworkControlWidget {
    constructor(
                network_figure,
                container_id
                )
    {
        var self = this;
        self.figure = network_figure;
                            
		var controlbox_width = 350,
			controlbox_height = 200,
			dy = 30,
            dx = 25;
		var
            n_grid_x = controlbox_width/dx,
			n_grid_y = controlbox_height/dy; // are used for putting a grid on the controls panels
	

		// this is the svg for the controls
		
	    self.controls = d3.select('#'+container_id).append("svg")
			.attr("width",controlbox_width)
			.attr("height",controlbox_height)
			.attr("class","explorable_widgets");
			//.style("border","1px solid black");
        //
        var defs = self.controls.append("defs");
        var gradient = defs.append("linearGradient");
        gradient.attr("id","animationSpeedGradient")
        /*
                .attr("x0","0")
                .attr("x1","1")
                .attr("y0","0")
                .attr("y1","0")
                */
            ;
        gradient.append("stop")
                .attr("offset","0%")
                .attr("stop-color", "#ff6464");
        gradient.append("stop")
               .attr("offset","50%")
                .attr("stop-color", "#fff");
        gradient.append("stop")
                .attr("offset","100%")
                .attr("stop-color", "#aaa");


		var g = widget.grid(controlbox_width,controlbox_height,n_grid_x,n_grid_y);
		var anchors = g.lattice(); // g has a method that returns a lattice with x,y coordinates

        /*
		// here we draw the lattice (usually not done in production)

		self.controls.selectAll(".grid").data(anchors).enter().append("circle")
			.attr("class","grid")
			.attr("transform",function(d){return "translate("+d.x+","+d.y+")"})
			.attr("r",1)
			.style("fill","black")
			.style("stroke","none");
        */
                        
		///////////////////
		// buttons
		///////////////////

		// we first define the button parameters

		var frame_left = { id: "frame_left", name:"←Frame", actions: ["rewind"], value: 0};
		var play_pause = { id: "play_pause", name:"", actions: ["play", "pause"], value: 1};
		//var record = { id: "record", name:"record", actions: ["record", "pause"], value: 1};
		var frame_right = { id: "frame_right", name:"Frame→", actions: ["reload"], value: 0};
		var stop = { id: "stop_and_reset", name:"", actions: ["stop"], value: 0};
		var capture = { id: "capture", name:"Pic", actions: ["capture"], value: 0};

        self.play_btn = widget.button(play_pause).size(60).symbolSize(30).update(function() { 
            if (self.figure.paused) 
            { 
                self.figure.play(); 
                self.hide_frame_buttons();
            } 
            else
            { 
                self.figure.pause();
                self.show_frame_buttons();
            } 
        });

        self.stop_btn = widget.button(stop).update(function() { 
                if (!self.figure.paused) 
                    self.play_btn.click();
                self.figure.global_it = 0;
                self.figure.update_all_visualizations();
            } ).size(40).symbolSize(20);

        self.frame_left_btn = widget.button(frame_left).size(40).symbolSize(20).update(function() { self.figure.update_global_it(-1); self.figure.update_all_visualizations(); } );
        self.frame_right_btn = widget.button(frame_right).size(40).symbolSize(20).update(function() { self.figure.update_global_it(+1); self.figure.update_all_visualizations(); } );
        self.capture_btn = widget.button(capture).size(40).symbolSize(20).update(function() { self.capture(); } );

		var buttons_main = [
            self.play_btn,
		];

		var buttonbox_main = g.block({x0:2, y0: n_grid_y - 2,width:1,height:0}).Nx(buttons_main.length);

		self.controls.selectAll(".buttons_main").data(buttons_main).enter().append(widget.buttonElement)
			.attr("transform",function(d,i){return "translate("+buttonbox_main.x(i)+","+buttonbox_main.y(0)+")"});	

		var buttons_frames = [
            self.stop_btn,
            self.frame_left_btn,
            self.frame_right_btn,
            self.capture_btn,
		];

		var buttonbox_frames = g.block({x0:4.5,y0: n_grid_y - 2,width:6,height:0}).Nx(buttons_frames.length);
		self.controls.selectAll(".buttons_frames").data(buttons_frames).enter().append(widget.buttonElement)
			.attr("transform",function(d,i){return "translate("+buttonbox_frames.x(i)+","+buttonbox_frames.y(0)+")"});	


        d3.select('body')
            .on("keydown", function() {
                var keyCode = d3.event.keyCode;
                if (keyCode == 37) // key: arrow left
                {
                    self.frame_left_btn.click();
                }
                else if (keyCode == 39) // key: arrow right
                {
                    self.frame_right_btn.click();
                }
                else if (keyCode == 32) // key: space
                {
                    self.play_btn.click();
                }
                else if (keyCode == 83 || keyCode == 124) // key: space
                {
                    self.capture_btn.click();
                }
            });
		///////////////////
		// sliders
		///////////////////	

		var speed_slider = { id: "play_speed", name: "Animation Speed", range: [-1,1], value: 0.6};
        console.log(self.figure.time);
		var aggregation_slider = { id: "aggregation_time", name: "Aggregation Time", range: [0.01,1], value: 0};


		var sliders = [
			widget.slider(speed_slider).label("bottom").update(function(){ 
                var val = speed_slider.value;
                var symbol, offset = +1;
                if (val < -0.6)
                    symbol = "◀◀" ;
                else if (val < 0)
                    symbol = "◀" ;
                else if (val < 0.6)
                    symbol = "▶" ;
                else
                    symbol = "▶▶" ;

                if (val < 0)
                    offset = -1;
                self.figure.update_animation_speed(speed_slider.value); 
                var handle_x = +d3.select("#slider_play_speed").select("circle").attr("cx");
                speed_label.attr("x",handle_x + offset)
                speed_label.attr("opacity",Math.sqrt(Math.abs(val)));
                speed_label.text(symbol);
                
            }),
			widget.slider(aggregation_slider).label("bottom").update(function(){ 
                self.figure.update_aggregation_time_scale(aggregation_slider.value);
            })
		]

		var sliderbox = g.block({x0:1,y0:n_grid_y-5.5,width:12.2,height:1.5}).Ny(sliders.length);

		sliders.forEach(function(d){
			d.width(sliderbox.w());
		})


		self.controls.selectAll(".slider").data(sliders).enter().append(widget.sliderElement)
			.attr("transform",function(d,i){return "translate("+sliderbox.x(0)+","+sliderbox.y(i)+")"});	

        ["#button_frame_right", "#button_frame_left"].forEach(function(id){
            d3.select(id)
                .style("opacity",0.3);
        });

        // add a gradient for
        d3.select("#slider_play_speed").select(".track-inset")
            .attr("y0",0)
            .attr("y1",0.0001) // this is needed because without
                               // changes in y, the coordinate transformation
                               // fails
            .style("stroke","url(#animationSpeedGradient)")
            .attr("class","none")
        ;

        // add a little play button

        var x_handle = +d3.select("#slider_play_speed").select("circle").attr("cx");
        var speed_label = d3.select("#slider_play_speed").insert("text",".track-overlay")
            .attr("class","tag")
            .attr("id","handle-label")
            .attr("text-anchor","tag")
            .attr("textLength","4")
            .attr("x",x_handle+1)
            .attr("y",4)
            .attr("text-anchor", "middle")
            .text("▶");

        sliders[0].click(0.6);
        //sliders[1].click(0);


   		///////////////////
		// toggles
		///////////////////
		// we first define the toggle parameters
		var record_toggle = {id:"recorder", name: "Record",  value: false};
		// now the array of toggle objets
		var toggles = [
			widget.toggle(record_toggle).label("bottom").update(function(d){
				var record_this = d.value;
				self.record(record_this);
			}),
		];

		// here comes the block for the toggles
		var togglebox = g.block({x0:12.2,y0:n_grid_y -2,width:1,height:1}).Ny(toggles.length);
		// and here we att them to the panel
		self.controls.selectAll(".toggle").data(toggles).enter().append(widget.toggleElement)
			.attr("transform",function(d,i){return "translate("+togglebox.x(i)+","+togglebox.y(i)+")"});	 
    }

    capture()
    {   
        var canvas = document.getElementById('mainCanvas');
        var dataURL = canvas.toDataURL("image/png");
        var max_num_length = (this.figure.time.length-1).toString().length;

        console.log(max_num_length);
        console.log(this.figure.global_it);
        d3.select("#gallery")
            .append("a")
                .attr("href",dataURL)
                .attr("target","_blank")
                .attr("download",
                        d3.format("0"+max_num_length+"d")
                                 (this.figure.global_it)
                    )
            .append("img")
                .attr("src",dataURL)
                .attr("width",150)
                .attr("class","screenshot");
    }

	record(start_recording)
	{
        var self = this;
        if(start_recording)
        {
            self.data = [];
            var canvas = document.querySelector("canvas")
            self.stream = canvas.captureStream(25);
            self.recorder = new MediaRecorder(self.stream, { mimeType: "video/webm" });
            self.recorder.ondataavailable = function(event) {
                if (event.data && event.data.size) {
                    self.data.push(event.data);
                }
            };
            self.recorder.onstop = () => {
                var url = URL.createObjectURL(new Blob(self.data, { type: "video/webm" }));
                self.downloadDataUrlFromJavascript("tacoma_recorded_video.webm",url);
            }
            self.recorder.start();
        }
        else
        {
            self.recorder.stop();
            delete self.recorder;
            delete self.stream;
        }

	}

    downloadDataUrlFromJavascript(filename, dataUrl) {

        // Construct the 'a' element
        var link = document.createElement("a");
        link.download = filename;
        link.target = "_blank";

        // Construct the URI
        link.href = dataUrl;
        document.body.appendChild(link);
        link.click();

        // Cleanup the DOM
        document.body.removeChild(link);
        //delete link;
    }

    hide_frame_buttons()
    {
        ["#button_frame_right", "#button_frame_left"].forEach(function(id){
            d3.select(id)
                .transition()
                .style("opacity",0.3);
        });
    }
    show_frame_buttons()
    {
        ["#button_frame_right", "#button_frame_left"].forEach(function(id){
            d3.select(id)
                .transition()
                .style("opacity",1);
        });
    }
}
