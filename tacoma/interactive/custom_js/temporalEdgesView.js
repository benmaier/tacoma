class temporalEdgesView {

    constructor(
                context,
                edges_data,
                plot_width,
                plot_height,
                offset_X,
                offset_Y,
                start_t,
                aggregation_time,
                fontsize_in_px,
                padding = 10,
                scale = 1,
                edge_line_width = 1,
                d3_format_string = "r"
                ) {


        this.context = context;
        this.data = edges_data;
        this.plot_width = plot_width;
        this.plot_height = plot_height;
        this.offset_X = offset_X;
        this.offset_Y = offset_Y;
        this.padding = padding;
        this.mark_T = start_t;
        this.edge_line_width = edge_line_width;
        this.fontsize = fontsize_in_px;
        this.scale = scale;
        this.d3_format_string = d3_format_string;
        

        this.aggregation_time = aggregation_time;

        this.xScale = d3.scaleLinear().range([this.offset_X + this.padding, 
                                              this.offset_X + this.plot_width-2*this.padding]);
        this.yScale = d3.scaleLinear().range([this.offset_Y + this.plot_height-2*this.padding, 
                                              this.offset_Y + this.padding]);

        this.edge_coordinates = [];

        this.critical_dx = this.plot_width / 5.0;
        this.critical_dy = this.plot_height / 3.0;
        this.fisheye = new fisheye(this.critical_dx,10,0.25);
        this.fisheye.range([this.offset_X + this.padding, this.offset_X + this.plot_width-this.padding]);

        // set domains of scales
        //
        this.t0 = this.data.xlim[0];
        this.tmax = this.data.xlim[1];
        this.t_half = d3.mean(this.data.xlim);

        this.xScale.domain( this.data.xlim );
        this.yScale.domain( this.data.ylim );

        var self = this;

        this.data.data.forEach( function (d) {
            var y = self.yScale(d[0]);
            var x1 = self.xScale(d[1]);
            var x2 = self.xScale(d[2]);
            self.edge_coordinates.push([ x1, x2, y ]);
        });

        this.draw(this.mark_T, -1, -1 );

        console.log(this.fisheye);

    }

    update_aggregation_time (aggregation_time)
    {
        this.aggregation_time = aggregation_time;
    }

    fisheye_X(x) {
        return this.fisheye.fisheye(x);
    }

    inverse_fisheye_X(x){
        return this.fisheye.fisheyeInverse(x);
    }
    /*

    fisheye_X(x, focusX) {
        var dx = x - focusX;
        if (Math.abs(dx) < this.critical_dx && Math.abs(dx) > 0)
        {
            var max_dx;
            if (dx > 0)
                max_dx = d3.min([this.plot_width - 2*this.padding - (focusX - this.offset_X), 
                                 this.critical_dx])
            else
                max_dx = d3.min([(focusX - this.offset_X) - this.padding, 
                                 this.critical_dx])

            return focusX + Math.sign(dx) * max_dx * Math.pow(Math.abs(dx)/max_dx,1/3);
        }
        else
            return x;
    }
    
    inverse_fisheye_X(X, focusX) {
        var dX = X - focusX;
        if (Math.abs(dX) < this.critical_dx && Math.abs(dX) > 0)
        {
            var max_dx;
            if (dX > 0)
                max_dx = d3.min([this.plot_width - 2*this.padding - (focusX - this.offset_X),
                                 this.critical_dx])
            else
                max_dx = d3.min([(focusX - this.offset_X) - this.padding,
                                 this.critical_dx])

            return focusX + Math.sign(dX) * max_dx * Math.pow(Math.abs(dX)/max_dx,3);
        }
        else
            return X;
    }

    fisheye_Y(y, mouseY) {
        var dy = y - mouseY;
        if (Math.abs(dy) < this.critical_dy)
        {
            var max_dy;
            if (dy > 0)
                max_dy = d3.min([this.plot_height - 2*this.padding - (mouseY - this.offset_Y), 
                                 this.critical_dy])
            else
                max_dy = d3.min([(mouseY - this.offset_Y) - this.padding, 
                                 this.critical_dy])

            return mouseY + Math.sign(dy) * max_dy * Math.sqrt(Math.abs(dy)/max_dy);
        }
        else
            return y;
    }
    */



    /*
        // this is a code block that deals with mouse events

        d3.select("canvas")
            .on('mouseout', function () {
                draw(this_mark_T,-1, -1);
            })
            .on('mousemove', function () {
                var mouseX = d3.event.layerX || d3.event.offsetX; 
                var mouseY = d3.event.layerY || d3.event.offsetY;
                draw(this_mark_T,mouseX, mouseY);
            })
            .on('click', function () {
                var mouseX = d3.event.layerX || d3.event.offsetX; 
                var mouseY = d3.event.layerY || d3.event.offsetY;

                if (mouseX > this.padding + this.offset_X && mouseX < this.plot_width - 2*this.padding + this.offset_X &&
                    mouseY > this.padding + this.offset_Y && mouseY < this.plot_height - 2*this.padding + this.offset_Y)
                {
                    this_mark_T = this.xScale.invert(mouseX);
                    draw(this.xScale.invert(mouseX), mouseX, mouseY);
                }
            });
    */

    draw(markT, mouseX, mouseY) {

        var mouse_in_picture = mouseX >= this.padding + this.offset_X && 
                      mouseY >= this.padding + this.offset_Y && 
                      mouseX <= this.plot_width - 2*this.padding + this.offset_X && 
                      mouseY <= this.plot_height - 2*this.padding + this.offset_Y;

        var fisheye = true;

        this.context.clearRect(this.offset_X,
                               this.offset_Y, 
                               this.plot_width, 
                               this.plot_height);
        this.context.fillStyle = "#fff";
        this.context.fillRect(this.offset_X,
                               this.offset_Y, 
                               this.plot_width, 
                               this.plot_height);

        this.context.fillStyle = 'rgb(0,0,0,1)';
        this.context.font = this.fontsize + "px Arial";

        var text;
        if (this.t0 == 0)
            text = "0";
        else
            text = d3.format(this.d3_format_string)(this.t0)
        this.context.textAlign = "left";
        this.context.fillText(text, this.offset_X + this.padding, 
                                    this.offset_Y + this.plot_height - this.fontsize/1.1);

        if (this.tmax == 0)
            text = "0";
        else
            text = d3.format(this.d3_format_string)(this.tmax)
        this.context.textAlign = "right";
        this.context.fillText(text, this.offset_X + this.plot_width - 2 * this.padding, 
                                    this.offset_Y + this.plot_height - this.fontsize/1.1);

        this.context.textAlign = "center";
        this.context.fillText("time", this.offset_X + this.plot_width / 2, 
                                      this.offset_Y + this.plot_height - this.fontsize/1.1);

        this.context.beginPath();
        this.context.strokeStyle = "rgb(100,100,100,0.8)";
        this.context.setLineDash([]);
        this.context.lineWidth = 1 * this.scale;
        this.context.moveTo(this.offset_X + this.plot_width - 2*this.padding, this.offset_Y + this.plot_height - 2*this.padding);
        this.context.lineTo(this.offset_X + this.padding, this.offset_Y + this.plot_height - 2*this.padding);
        this.context.lineTo(this.offset_X + this.padding, this.offset_Y + this.padding);
        this.context.stroke();
        this.context.setLineDash([]);

        var markX;
        if (markT>=0)
        {
            markX = this.xScale(markT);
            //mouseX = markX;
            this.fisheye.focus(markX);
        }

        this.context.strokeStyle = 'rgba(0, 0, 0, 1)';
        this.context.lineWidth = this.edge_line_width;
        this.context.beginPath();


        for(var i = 0;  i < this.edge_coordinates.length; i++)
        {
            var d = this.edge_coordinates[i];
            var x1 = d[0];
            var x2 = d[1];
            var y = d[2];

            if (fisheye)
            {
                var old_x1 = x1;
                x1 = this.fisheye_X(x1);
                x2 = this.fisheye_X(x2);
                //x1 = 0.5*(this.fisheye_X(x1, mouseX) + this.fisheye_X(x1, markX));
                //x2 = 0.5*(this.fisheye_X(x2, mouseX) + this.fisheye_X(x2, markX));
                //x2 = this.fisheye_X(x2, mouseX);
                //y = this.fisheye_Y(y, mouseY);
            }

            this.context.moveTo(x1,y);
            this.context.lineTo(x2,y);
        }
        this.context.closePath();
        this.context.stroke();

        if (markT >= 0)
        {
            this.context.beginPath();
            //var markX = this.xScale(markT);
            var markX2 = this.xScale(markT+this.aggregation_time);
            if (fisheye)
            {
                markX = this.fisheye_X(markX);
                markX2 = this.fisheye_X(markX2);
                //markX = this.fisheye_X(markX, mouseX);
                //markX2 = this.fisheye_X(markX2, mouseX);
            }
            markX2 = d3.min([markX2, this.offset_X + this.plot_width - 2*this.padding ]);
            var rect_width = markX2 - markX;
            var rect_height = Math.floor(this.plot_height - 2*this.padding);
            var rect_y = Math.floor(this.offset_Y)+0.5*this.padding;

            this.context.lineWidth = 1.5*this.scale;
            this.context.strokeStyle = 'rgb(150,150,150,0.8)';
            this.context.rect(markX, rect_y, rect_width, rect_height);
            this.context.closePath();
            this.context.fillStyle = 'rgb(150,150,150,0.2)';
            this.context.fill();
            this.context.stroke();

            this.context.beginPath();
            this.context.lineWidth = 2*this.scale;
            //this.context.setLineDash([10,20]);
            this.context.strokeStyle = 'rgb(255,100,100,0.8)';
            this.context.moveTo(markX, this.offset_Y);
            this.context.lineTo(markX, this.plot_height + this.offset_Y - this.padding);
            this.context.closePath();
            this.context.stroke();

        }

        var markX_new = null;
        if (mouse_in_picture)
        {
            markX_new = this.inverse_fisheye_X(mouseX);

            //var markX = this.xScale(markT);
            var this_time = this.xScale.invert(markX_new);
            var markX2 = this.xScale(this_time+this.aggregation_time);
            if (fisheye)
            {
                markX_new = this.fisheye_X(markX_new);
                markX2 = this.fisheye_X(markX2);
                //markX = this.fisheye_X(markX, mouseX);
                //markX2 = this.fisheye_X(markX2, mouseX);
            }
            markX2 = d3.min([markX2, this.offset_X + this.plot_width - 2*this.padding ]);
            var rect_width = markX2 - markX_new;
            var rect_height = Math.floor(this.plot_height - 2*this.padding);
            var rect_y = Math.floor(this.offset_Y)+0.5*this.padding;

            this.context.beginPath();
            this.context.lineWidth = 1.5*this.scale;
            this.context.setLineDash([3*this.scale,3*this.scale]);
            this.context.strokeStyle = 'rgb(150,150,150,0.8)';
            this.context.rect(markX_new, rect_y, rect_width, rect_height);
            this.context.fillStyle = 'rgb(150,150,150,0.2)';
            this.context.stroke();
            this.context.fill();


            this.context.beginPath();
            this.context.lineWidth = 3*this.scale;
            //this.context.setLineDash([10,20]);
            this.context.setLineDash([10*this.scale,20*this.scale]);
            this.context.strokeStyle = 'rgb(255,100,100,0.8)';
            this.context.moveTo(markX_new, this.offset_Y);
            this.context.lineTo(markX_new, this.plot_height + this.offset_Y - this.padding);
            this.context.closePath();
            this.context.stroke();

            this.context.save();
            this.context.fillStyle = 'rgb(255,100,100,0.8)';
            //this.context.translate(markX_new, this.plot_height / 2 + this.offset_Y);
            //this.context.rotate(-Math.PI/2);
            this.context.font = this.fontsize + "px Arial";
            var text = " " + d3.format(this.d3_format_string)(this_time) + " ";
            var dx = markX - markX_new;
            if (dx > 0 && dx < text.length*this.fontsize)
                this.context.textAlign = 'right';
            else
                this.context.textAlign = 'left';
            this.context.fillText(text, markX_new, this.offset_Y + this.fontsize / 1.3);
            this.context.restore()
        }

        var text = "   " + d3.format(this.d3_format_string)(this.xScale.invert(markX)) + " ";
        if (!(markX_new === null))
        {
            var dx = markX - markX_new;
            if (dx > 0 && dx < text.length*this.fontsize)
                this.context.textAlign = 'left';
            else
                this.context.textAlign = 'right';
        }
        else if (markX - this.offset_X < text.length*this.fontsize / 2)
            this.context.textAlign = 'left';
        else
            this.context.textAlign = 'right';

        this.context.fillStyle = 'rgb(255,100,100,0.8)';
        //this.context.translate(markX, this.plot_height / 2 + this.offset_Y);
        //this.context.rotate(-Math.PI/2);
        this.context.font = this.fontsize + "px Arial";
        this.context.fillText(text, markX, this.offset_Y + this.fontsize/1.3);

        this.context.save();
        this.context.fillStyle = 'rgb(0,0,0,1)';
        this.context.translate(this.offset_X, this.offset_Y + (this.plot_height-this.padding) / 2);
        this.context.rotate(-Math.PI/2);
        this.context.textAlign = "center";
        this.context.fillText("active edges", 0, this.fontsize/1.5);
        this.context.restore();

    }
}
