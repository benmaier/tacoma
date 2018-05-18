var widget = (function(){

	var widget; // instance of widget
	
	// module stuff
	
	var grid = function(W,H,Nx,Ny){
		
	  	var X = d3.scaleLinear().domain([0,Nx]).range([0,W]);
		var Y = d3.scaleLinear().domain([0,Ny]).range([H,0]);
		
		var lattice = function(){
			return d3.range((Nx+1)*(Ny+1)).map(function(i){
   				return { m:(i % (Nx+1)), n: Math.floor(i / (Nx+1)), x: X((i % (Nx+1))), y: Y(Math.floor(i / (Nx+1)))}
   	 		})		
		}
		

		var block = function(position){
			
			var x0,y0,w,h,nx,ny;
			var edge = "[]";		
			
			if ("undefined" === typeof position) {
				x0 = 0; y0 = 0; w = 1; h = 1;
			} else {
				x0 = position.x0;
				y0 = position.y0;
				w = position.width;
				h = position.height;				
			}
			
			nx = w+1; ny = h+1;
			
			var x = function(k){ return nx == 1 ? X(x0) : X(x0 + k * w/(nx-1)) };
			var y = function(k){ return ny == 1 ? Y(y0) : Y(y0 + k * h/(ny-1)) };
			
			var lattice = function(){
				return d3.range((nx)*(ny)).map(function(i){
	   				return { 
						m:(i % (nx)), 
						n: Math.floor(i / (nx)), 
						x: x((i % (nx))), 
						y: y(Math.floor(i / (nx)))}
	   	 		})		
			}
	
			// this is what block returns
			return {				
				x0 : function(arg) { if ("undefined" === typeof arg) { return x0 } else { x0 = arg; return this }}, 
				y0 : function(arg) { if ("undefined" === typeof arg) { return y0 } else { y0 = arg; return this }},
				width : function(arg) { if ("undefined" === typeof arg) { return w } else { w = arg ; return this}},
				height : function(arg) { if ("undefined" === typeof arg) { return h } else { h = arg ; return this}},	
				Nx : function(arg) { if ("undefined" === typeof arg) { return nx } else { nx = arg ; return this}},
				Ny : function(arg) { if ("undefined" === typeof arg) { return ny } else { ny = arg ; return this}},					
				x:x, 
				y:y, 
				w:function(){return Math.abs(X(w+x0)-X(x0))},
				h:function(){return Math.abs(Y(h+y0)-Y(y0))},
				position: function(){return {x0:x0,y0:y0,width:w,height:h}},
				lattice:lattice
			}
		}
		
		// this is what grid returns
		
		return { x:X, y:Y, block:block, lattice:lattice }
							
	};	
	var button = function(parameter){

		var size = 50,
			symbolSize = 25,
			shape = "round",
			label = "bottom",
			fontSize = 12,

			update = function(x){};
			click = function (){
				
				parameter.value = (parameter.value + 1) % parameter.actions.length;
			
				d3.select("#button_"+parameter.id).select("class",".button-background")
					.transition().duration(1000).attr("class","button-background-lit")
					.transition().duration(1000).attr("class","button-background")
			
				d3.select("#button_"+parameter.id).selectAll("path")
					.attr("d",symbol(parameter.actions[parameter.value],symbolSize/2))
					.transition().attr("class","button-symbol-lit")
					.transition().attr("class","button-symbol")
			
				update(parameter);
			};
			
			return {
				size: function(arg) { if ("undefined" === typeof arg) { return size } else { size = arg; return this }},
				symbolSize: function(arg) { if ("undefined" === typeof arg) { return symbolSize } else { symbolSize = arg; return this }},
				shape: function(arg) { if ("undefined" === typeof arg) { return shape } else { shape = arg; return this }},
				label: function(arg) { if ("undefined" === typeof arg) { return label } else { label = arg; return this }},
				fontSize: function(arg) { if ("undefined" === typeof arg) { return fontSize } else { fontSize = arg; return this }},
				parameter: function(arg) { if ("undefined" === typeof arg) { return parameter } else { parameter = arg; return this }},
				name: function() {return  parameter.name},
				id: function() {return parameter.id},
				value: function() {return parameter.value},
				actions : function() {return parameter.actions},
				update: function(arg) { if ("function" === typeof arg) {update = arg; return this} else { update(arg) }},
				click:click
			}
	};	
	var buttonElement = function(d,i){

		var backbox ;
		var hakenschniepel = document.createElementNS("http://www.w3.org/2000/svg", "g");	

		var s = d3.select(hakenschniepel).append("g")
			.attr("class", "button")
			.attr("id", "button_" + d.id())

		if (d.shape()=="rect"){
			backbox = s.append("rect")
				.attr("width",d.size())
				.attr("height",d.size())
				.attr("transform","translate("+(-d.size()/2)+","+(-d.size()/2)+")")
				.attr("rx",5).attr("ry",5)
		} else {
			backbox = s.append("circle")
				.attr("r",d.size()/2)	
		}
	
		backbox.attr("class","button-background")
			.on("mouseover",function(x){
				d3.select(this).attr("class","button-background-hover")
				d3.select("#button_" + d.id()).select("path")
					.attr("class","button-symbol-hover")
			})
			.on("mouseout",function(){
				d3.select(this).attr("class","button-background")
				d3.select("#button_" + d.id()).select("path")
					.attr("class","button-symbol")
			})
			.on("click",d.click)
	
		s.append("path")
			.attr("d",symbol(d.actions()[d.value()],d.symbolSize()/2))
			.attr("class","button-symbol")
	
	
		if (d.name!=""){
		
				var xpos,ypos,anchor,valign;
		
				if(d.label() == "top") {
					xpos = 0;
					ypos = -d.size()/2-5;
					anchor = "middle";
					valign ="bottom";
			
				}
		
				if(d.label() == "bottom") {
					xpos = 0;
					ypos = d.size()/2+5;
					anchor = "middle";
					valign = "hanging";
				}
		
				if(d.label() == "right") {
					xpos = d.size() / 2 + 5;
					ypos = 0;
					anchor = "start";
					valign ="middle";
				}
		
				if(d.label() == "left") {
					xpos =  - d.size() / 2 - 5;
					ypos = 0;
					anchor = "end";
					valign ="middle";
				}

				s.append("text").text(d.name())
					.attr("class", "tag")
					.style("opacity",1)
					.style("text-anchor",anchor)
					.style("font-size",d.fontSize())
					.style("alignment-baseline",valign)
					.attr("transform", "translate(" + (xpos) + "," + (ypos) + ")")
		
			}
	
	
		return hakenschniepel;
	};	
	var toggle = function(parameter){
		var size = 9,
			fontSize = 12,
			border = 0.5,
			label = "top",
			update = function(x){};

			var click = function(){
				parameter.value = ! parameter.value;
				d3.selectAll("#handle_" + parameter.id).transition()
					.attr("cx", parameter.value ? 2*size : 0)
				d3.selectAll("#inset_"+parameter.id)
					.attr("class",parameter.value ?  "track-inset-lit" : "track-inset")
				update(parameter);
			}
			
			return {
				size: function(arg) { if ("undefined" === typeof arg) { return size } else { size = arg; return this }},
				border: function(arg) { if ("undefined" === typeof arg) { return border } else { border = arg; return this }},
				label: function(arg) { if ("undefined" === typeof arg) { return label } else { label = arg; return this }},
				fontSize: function(arg) { if ("undefined" === typeof arg) { return fontSize } else { fontSize = arg; return this }},
				parameter: function(arg) { if ("undefined" === typeof arg) { return parameter } else { parameter = arg; return this }},
				name: function() {return  parameter.name},
				id: function() {return parameter.id},
				value: function() {return parameter.value},
				update: function(arg) { if ("function" === typeof arg) {update = arg; return this} else { update(arg) }},
				click:click
			}
			
	}	
	var toggleElement = function(d,i){
	
	
		var hakenschniepel = document.createElementNS("http://www.w3.org/2000/svg", "g");	
	
		var s = d3.select(hakenschniepel)
			.attr("class", "toggle")
			.attr("id", "toggle_" + d.id())

		s.append("line")
			.attr("class", "track")
			.attr("x1", 0).attr("x2", 2*d.size())
			.style("stroke-width", 2*d.size() + 2 * d.border())
	
		s.append("line")
			.attr("id", "inset_" + d.id())
			.attr("class", d.parameter().value ? "track-inset-lit" : "track-inset")
			.attr("x1", 0).attr("x2", 2*d.size())
			.style("stroke-width", 2*d.size())
	
		s.append("line")
			.attr("class", "track-overlay")
			.attr("x1", 0).attr("x2", 2*d.size())
			.style("stroke-width", 2* d.size())
			.on("click",d.click)
		
		s.insert("circle", ".track-overlay")
			.attr("class", "handle")
			.attr("id", "handle_" + d.id())
			.attr("r", d.size())
			.attr("cx", d.value() ? 2*d.size() : 0);

		if (d.name()!=""){
		
			var xpos,ypos,anchor,valign;
		
			if(d.label() == "top") {
				xpos = d.size();
				ypos = -2*d.size();
				anchor = "middle";
				valign ="middle";
			
			}
		
			if(d.label() == "bottom") {
				xpos = d.size();
				ypos = 2*d.size();
				anchor = "middle";
				valign ="middle";
			}
		
			if(d.label() == "right") {
				xpos = d.size() + 3 * d.size();
				ypos = 0;
				anchor = "start";
				valign ="middle";
			}
		
			if(d.label() == "left") {
				xpos = d.size() -3 * d.size();
				ypos = 0;
				anchor = "end";
				valign ="middle";
			}

			s.append("text").text(d.name())
				.attr("class", "tag")
				.style("opacity",1)
				.style("text-anchor",anchor)
				.style("font-size",d.fontSize())
				.style("alignment-baseline",valign)
				.attr("transform", "translate(" + (xpos) + "," + (ypos) + ")")
		
		}
	

		return hakenschniepel;	
	
	
	}	
	var slider = function(parameter){
	
		var width = 100,
			handleSize = 10,
			trackSize = 6,
			fontSize = 12,
			trackBorder = 0.5,
			label = "top-left", 
			update = function(x){};
			
		
	
		var click = function(x) {
			
			var X = d3.scaleLinear().domain(parameter.range).range([0, width]).clamp(true);
						d3.selectAll("#handle_" + parameter.id).transition().attr("cx", X(x))
						parameter.value = x;
						update();
		}
		
		
			
		return {
			width: function(arg) { if ("undefined" === typeof arg) { return width } else { width = arg; return this }},
			handleSize: function(arg) { if ("undefined" === typeof arg) { return handleSize } else { handleSize = arg; return this }},
			trackSize: function(arg) { if ("undefined" === typeof arg) { return trackSize } else { trackSize = arg; return this }},
			trackBorder: function(arg) { if ("undefined" === typeof arg) { return trackBorder } else { trackBorder = arg; return this }},
			label: function(arg) { if ("undefined" === typeof arg) { return label } else { label = arg; return this }},
			fontSize: function(arg) { if ("undefined" === typeof arg) { return fontSize } else { fontSize = arg; return this }},
			parameter: function(arg) { if ("undefined" === typeof arg) { return parameter } else { parameter = arg; return this }},
			name: function() {return  parameter.name},
			id: function() {return parameter.id},
			range: function() {return parameter.range},
			value: function() {return parameter.value},
			update: function(arg) { if ("function" === typeof arg) {update = arg; return this} else { update(arg) }},
			click:click
		}
		
	}	
	var sliderElement = function(d,i) {
	
		var hakenschniepel = document.createElementNS("http://www.w3.org/2000/svg", "g")

		d.X = d3.scaleLinear()
			.domain(d.range())
			.range([0, d.width()]).clamp(true);

	
		var s = d3.select(hakenschniepel)
			.attr("class", "slider")
			.attr("id", "slider_" + d.id())

		s.append("line")
			.attr("class", "track")
			.attr("x1", d.X.range()[0]).attr("x2", d.X.range()[1])
			.style("stroke-width", d.trackSize() + 2 * d.trackBorder())

		s.append("line")
			.attr("class", "track-inset")
			.attr("x1", d.X.range()[0]).attr("x2", d.X.range()[1])
			.style("stroke-width", d.trackSize())

		s.append("line")
			.attr("class", "track-overlay")
			.attr("x1", d.X.range()[0]).attr("x2", d.X.range()[1])
			.style("stroke-width", 2* d.handleSize())
			
			.call(d3.drag()
				.on("start drag", function() {
					var value = d.X.invert(d3.event.x);
					d3.selectAll("#handle_" + d.id()).attr("cx", d.X(value))
					d.parameter().value = value;
					d.update(d);
				})
				
			);

		s.insert("circle", ".track-overlay")
			.attr("class", "handle")
			.attr("id", "handle_" + d.id())
			.attr("r", d.handleSize())
			.attr("cx", d.X(d.value()));


		if(d.name()!=""){					
	
			var xpos,ypos,anchor,valign="middle";
		
			ypos = d.label().match(/bottom/i)!=null ? 2 * d.handleSize() : - 2 * d.handleSize();
			xpos = d.label().match(/right/i)!=null ? d.width() : (d.label().match(/center/i)!=null ? d.width() / 2 : 0);
			anchor = d.label().match(/right/i)!=null ? "end" : (d.label().match(/center/i)!=null ? "middle" : "start")
				
			s.append("text").text(d.name())
				.attr("class", "tag")
				.style("text-anchor",anchor)
				.style("alignment-baseline",valign)
				.style("font-size",d.fontSize())
				.style("opacity",1)
				.attr("transform", "translate(" + (xpos) + "," + (ypos) + ")")
		}
	
		return hakenschniepel;
	}	
	var radio = function (parameter){
	var size = 200,
		buttonSize = 20,
		buttonInnerSize = 12,
		fontSize = 12,
		padding = 5,
		orientation = "vertical",
		shape = "rect",
		label = "right",
		update = function(x){};
	
		var click = function(j){
				parameter.value=j;
				d3.select("#radio_"+parameter.id).selectAll(".radiobutton-on").attr("class","radiobutton-off");
				d3.select("#radio_"+parameter.id).selectAll(".radiobutton-off").attr("class",function(z,k){
					return k==j ? "radiobutton-on" : "radiobutton-off"
				});
				update();
		}	
		
		return {
			size: function(arg) { if ("undefined" === typeof arg) { return size } else { size = arg; return this }},
			buttonSize: function(arg) { if ("undefined" === typeof arg) { return buttonSize } else { buttonSize = arg; return this }},
			buttonInnerSize: function(arg) { if ("undefined" === typeof arg) { return buttonInnerSize } else { buttonInnerSize = arg; return this }},
			padding: function(arg) { if ("undefined" === typeof arg) { return padding} else { padding = arg; return this }},
			orientation: function(arg) { if ("undefined" === typeof arg) { return orientation} else { orientation = arg; return this }},
			shape: function(arg) { if ("undefined" === typeof arg) { return shape } else { shape = arg; return this }},
			label: function(arg) { if ("undefined" === typeof arg) { return label } else { label = arg; return this }},
			fontSize: function(arg) { if ("undefined" === typeof arg) { return fontSize } else { fontSize = arg; return this }},
			parameter: function(arg) { if ("undefined" === typeof arg) { return parameter } else { parameter = arg; return this }},
			name: function() {return  parameter.name},
			id: function() {return parameter.id},
			value: function() {return parameter.value},
			choices : function() {return parameter.choices},
			update: function(arg) { if ("function" === typeof arg) {update = arg; return this} else { update(arg) }},
			click:click
		}
	}
	var radioElement = function(d,i){
	
		var hakenschniepel = document.createElementNS("http://www.w3.org/2000/svg", "g");
		var N = d.choices().length;
		var n = d3.range(N);
		var X = d3.scaleLinear().domain([0,N]).range([-d.size(),0]);
	
		var checkbox = d3.select(hakenschniepel).attr("class","radio").attr("id","radio_"+d.id());

		var button = checkbox.selectAll(".radiobutton").data(n).enter().append("g")
			.attr("class","radiobutton")
			.attr("id",function(x,j){return "radiobutton_"+ d.id() + "_" +j})
			.attr("transform",function(x,j){
				return d.orientation()=="vertical" ? "translate(0,"+X(j)+")" : "translate("+X(j)+",0)";
			})
	
		var background, led;
	
		if (d.shape()=="rect"){
		
			background = button.append("rect")
				.attr("width",d.buttonSize())
				.attr("height",d.buttonSize())
				.attr("rx",2)
				.attr("ry",2)
				.attr("class","radiobutton-background")
				.attr("transform","translate("+(-d.buttonSize()/2)+","+(-d.buttonSize()/2)+")")
		
			led = button.append("rect")
				.attr("width",d.buttonInnerSize())
				.attr("height",d.buttonInnerSize())
				.attr("rx",2)
				.attr("ry",2)
				.attr("class","radiobutton-off")
				.attr("transform","translate("+(-d.buttonInnerSize()/2)+","+(-d.buttonInnerSize()/2)+")")
				.attr("class",function(x,j){return j==d.value() ? "radiobutton-on" : "radiobutton-off"})
		} else {
		
			background = button.append("circle")
				.attr("r",d.buttonSize()/2)	
				.attr("class","radiobutton-background")
		
			led = button.append("circle")
				.attr("r",d.buttonInnerSize()/2)
				.attr("class",function(x,j){return j==d.value() ? "radiobutton-on" : "radiobutton-off"})
		
		}
	
		background
			.on("mouseover",function(x,j){
				d3.select("#radiobutton_"+d.id() + "_" +j).select(".radiobutton-off")
				.attr("class","radiobutton-hover")

			})
			.on("mouseout",function(){
					led.attr("class",function(x,j){return j==d.value() ? "radiobutton-on" : "radiobutton-off"})
			})
			.on("click",function(x,j){
				d.parameter().value=j;
				led.attr("class",function(z,k){return k==d.value() ? "radiobutton-on" : "radiobutton-off"})
				d.update(d);
			})
	
	
		var xpos,ypos,anchor,valign;

		xpos = d.buttonSize();
		ypos = 0;
		anchor = "start";
		valign = "middle";
	
		if(d.label()=="left"){
			xpos = -d.buttonSize();
			ypos = 0;
			anchor = "end";
			valign = "middle";
		}
	
		if(d.label()=="bottom"){
			ypos = d.buttonSize();
			xpos = 0;
			anchor = "middle";
			valign = "hanging";
		}
	
		if(d.label()=="top"){
			ypos = -d.buttonSize();
			xpos = 0;
			anchor = "middle";
			valign = "bottom";
		}
	
	
		button.append("text")
			.attr("class","tag")
			.text(function(x,j){return d.choices()[j]})
			.attr("alignment-baseline",valign)
			.attr("transform","translate("+(xpos)+","+ypos+")")
			.style("font-size",d.fontSize())
			.attr("text-anchor",anchor)
	
		return hakenschniepel;
	}
	
	
	return {
			grid: function(W,H,M,N){
				widget = grid(W,H,M,N);
				return widget;
			},
			button: function(arg){
				return button(arg);
			},
			toggle: function(arg){
				return toggle(arg);
			},
			slider: function(arg){
				return slider(arg);
			},
			radio: function(arg){
				return radio(arg);
			},
			buttonElement: function(d,i){
				return buttonElement(d,i);
			},
			toggleElement: function(d,i){
				return toggleElement(d,i);
			},
			sliderElement: function(d,i){
				return sliderElement(d,i);
			},
			radioElement: function(d,i){
				return radioElement(d,i);
			}
		}	
	
	function symbol(type,scale){
	
			if (typeof scale === "undefined") {scale = 100}
	
			switch (type) {
			case "play":
				return function() {
						var p = d3.path();
						p.moveTo(scale * 1, scale * 0);
						p.lineTo(scale * (-0.5), scale * (Math.sqrt(3) / 2))
						p.lineTo(scale * (-0.5), scale * ( - Math.sqrt(3) / 2))
						p.closePath();
	
						return p.toString();
					}	
					break;
			case "back":
				return function() {
						var p = d3.path();
						p.moveTo( - scale * 1, scale * 0);
						p.lineTo(scale * (0.5), scale * (Math.sqrt(3) / 2))
						p.lineTo(scale * (0.5), scale * ( - Math.sqrt(3) / 2))
						p.closePath();

						return p.toString();
					}	
						break;		
			case "pause":
					return function() {
							var g = 1 / 3;
							var p = d3.path();
							var c = 0.9
							p.moveTo(scale * c, scale * c);
							p.lineTo(scale * c, scale * (-c))
							p.lineTo(scale * (c * g), scale * ( - c))
							p.lineTo(scale * (c * g), scale * (  c))
							p.closePath();
		
							p.moveTo(- scale * c, scale * c);
							p.lineTo(- scale * c, scale * (-c))
							p.lineTo(- scale * (c * g), scale * ( - c))
							p.lineTo(- scale * (c * g), scale * (  c))
							p.closePath();
		
	
							return p.toString();
						};
						break;
				case "reload":
					return function() {
		
							var theta = Math.PI/2.5;
							var theta1 = theta / 2;
							var theta0 = 2*Math.PI - theta / 2;
							var width = 0.5;
							var arrow_width = 0.6;
							var arrow_height = 0.6;
		
							var p = d3.path();
		
							p.moveTo(scale * Math.cos (theta0), scale * Math.sin(theta0));
							p.arc(0,0,scale,theta0,theta1,true);
							p.lineTo(scale *(1-width) * Math.cos (theta1), scale *(1-width) * Math.sin (theta1))
							p.arc(0,0,scale * (1-width),theta1,theta0,false);
							p.lineTo(scale * (1 - arrow_width - width / 2 ) * Math.cos(theta0), scale * (1 - arrow_width - width / 2  ) * Math.sin(theta0))
		
							var w0 = [scale *(1 - width / 2) * Math.cos(theta0),scale * (1 - width / 2) * Math.sin(theta0)]
							var z = [scale * arrow_height * Math.cos(theta0+Math.PI / 2), scale * arrow_height * Math.sin(theta0+Math.PI / 2)] 
		
							p.lineTo(w0[0]+z[0], w0[1]+z[1])
							p.lineTo(scale * (1 + arrow_width - width / 2  ) * Math.cos(theta0), scale * (1 + arrow_width - width / 2 ) * Math.sin(theta0))
		
		
		
							p.closePath();
		
	
							return p.toString();
						};
						break;
                    case "record":
                        return function() {
                            var p = d3.path();
                            p.arc(0,0,scale,0,2*Math.PI);
                            p.closePath();
                            return p;
                        }
                        break;
                    case "capture":
                        return function() {
                            var p = d3.path();
                            var theta0 = Math.PI/10;
                            var sign = +1;
                            var inner_radius = scale / 2;
                            var angle1 = Math.PI - theta0;
                            var angle2 = theta0;
                            var angle3 = -theta0;
                            var angle4 = Math.PI + theta0;

                            p.arc(0,0,inner_radius, angle4, angle3);                            
                            p.lineTo(scale, inner_radius*Math.sin(angle4));
                            p.lineTo(scale, -scale);
                            p.lineTo(-scale, -scale);
                            p.lineTo(-scale, inner_radius*Math.sin(angle4));

                            p.closePath();

                            p.arc(0,0,inner_radius, angle2, angle1);
                            p.lineTo(-scale, inner_radius*Math.sin(angle1));
                            p.lineTo(-scale, scale);
                            p.lineTo(scale, scale);
                            p.lineTo(scale, inner_radius*Math.sin(angle1));

                            p.closePath();

                            return p;
                        }
                        break;
					case "rewind":
						return function() {
		
								var theta = Math.PI/2.5;
								var theta0 = theta / 2+ Math.PI;
								var theta1 = 2*Math.PI - theta / 2 + Math.PI;
								var width = 0.5;
								var arrow_width = 0.6;
								var arrow_height = -0.6;
		
								var p = d3.path();
		
								p.moveTo(scale * Math.cos (theta0), scale * Math.sin(theta0));
								p.arc(0,0,scale,theta0,theta1,false);
								p.lineTo(scale *(1-width) * Math.cos (theta1), scale *(1-width) * Math.sin (theta1))
								p.arc(0,0,scale * (1-width),theta1,theta0,true);
								p.lineTo(scale * (1 - arrow_width - width / 2 ) * Math.cos(theta0), scale * (1 - arrow_width - width / 2  ) * Math.sin(theta0))
		
								var w0 = [scale *(1 - width / 2) * Math.cos(theta0),scale * (1 - width / 2) * Math.sin(theta0)]
								var z = [scale * arrow_height * Math.cos(theta0+Math.PI / 2), scale * arrow_height * Math.sin(theta0+Math.PI / 2)] 
		
								p.lineTo(w0[0]+z[0], w0[1]+z[1])
								p.lineTo(scale * (1 + arrow_width - width / 2  ) * Math.cos(theta0), scale * (1 + arrow_width - width / 2 ) * Math.sin(theta0))
		
		
		
								p.closePath();
		
	
								return p.toString();
							};
							break;		
						
				case "stop":
						return 	function() {
				var p = d3.path();
				var c = 0.9
				p.moveTo(scale * c, scale * c);
				p.lineTo(scale *(-c), scale * c)
				p.lineTo(scale * (-c), scale * ( - c))
				p.lineTo(scale * (c), scale * ( - c))
				p.closePath();
	
				return p.toString();
			}	;
			break;
			default:
				return function(){
					var p = d3.path();
					p.arc(0,0,scale,0,2*Math.PI,true);
					p.closePath();
					return p.toString();
				};
		
					
			}
		}
	
	
})()
