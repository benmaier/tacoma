class fisheye {
    constructor (
        r,
        d = 3,
        xw = 0.4,
        mode = 'continuous'
    )
    {
        this.r = r;
        this.x_max = null;
        this.x_min = null;

        this.range_x = null;
        this.range_y = null;
        this.magnification(d)
            .demagnificationWidth(xw);
        this.mode = mode;

        this.focus_x = null;
    }

    fisheye(x,direction = 0, inverse = false){

        if (this.focus_x === null)
            return x;
      
        if ((this.xw == 1 || this.d == 0) && this.mode == 'continuous')
            return x;
      
        let focus;
        if (!this.focus_x.length)
          focus = this.focus_x;
        else
          focus = this.focus_x[direction];

        let dx = x - focus;
      
        if ((Math.abs(dx) > this.r) || (x == focus))
            return x;


        let dmax = this.r;
        let range;
        if (direction == 0)
            range = this.range_x;
        else
            range = this.range_y;

        if (dx > 0){
            if (! (range === null) )
                if ( range[1] - focus < this.r)
                    dmax = range[1] - focus;
        } 
        else if (dx < 0)
        {
            if (! (range === null) )
                if ( focus - range[0] < this.r)
                    dmax = focus - range[0];
        }

        let d = this.d;
        let xc = this.xc;
        let A1 = this.A1;
        let A2 = this.A2;

        let rescaled = Math.abs(dx) / dmax;

        let new_r;

        if (this.mode == 'continuous')
        {
            if (!inverse)
              new_r = this.fisheyeContinuous(rescaled);
            else
              new_r = this.fisheyeContinuousInverse(rescaled);
        }
        else if (this.mode == 'sarkar')
        {
            if (!inverse)
              new_r = this.fisheyeSarkar(rescaled);
            else
              new_r = this.fisheyeSarkarInverse(rescaled);
        }
        else
            throw 'Unknown mode "'+ this.mode + '"';

        return focus + Math.sign(dx) * dmax * new_r;
    }

    fisheyeInverse(x,direction=0){
        return this.fisheye(x,direction,true);
    }
  
    fisheyeFunction(x){
      if ( (x<=0) || (x>=1) )
        return x;
      
      let f;
      if (this.mode == 'continuous')
        f = this.fisheyeContinuous(x);
      else if (this.mode == 'sarkar')
        f = this.fisheyeSarkar(x);     
      
      return f;
    }
  
    fisheyeInverseFunction(x){
      if ( (x<=0) || (x>=1) )
        return x;
      
      let f;
      if (this.mode == 'continuous')
        f = this.fisheyeContinuousInverse(x);
      else if (this.mode == 'sarkar')
        f = this.fisheyeSarkarInverse(x);     
      
      return f;
    }

    fisheyeCartesian(pos,inverse=false) {
        let x = this.fisheye(pos[0],0,inverse);
        let y = this.fisheye(pos[1],1,inverse);
        return [x,y];
    }

    fisheyeCartesianInverse(x) {
        return this.fisheyeCartesian(x, true);
    }

    solve_linear_equation(a,b,c,d,vector_b) {
        /* solves equation 
          / a  b \  / A1 \   / b1 \
          |      |  |    | = |    |
          \ c  d /  \ A2 /   \ b2 /
        */
        let det = a*d - b*c;
        let A1 = ( d*vector_b[0] - b*vector_b[1]) / det;
        let A2 = (-c*vector_b[0] + a*vector_b[1]) / det;
        return [ A1, A2 ];
    }

    rescale()
    {
        let xw = this.xw;
        let d = this.d;
        if ((xw>0) && (xw<1))
        {
          let A = this.solve_linear_equation( Math.pow(xw,2)/2., 1 - ((d+1)*xw / (d*xw+1)), 
                                              xw,                - (d+1) / Math.pow(d*xw+1,2),                                                                   [ (1-xw), -1 ]);
                   
          this.A1 = A[0];
          this.A2 = A[1];
        } 
        else if (xw == 0)
        {
          this.A1 = 0;
          this.A2 = 1;
        }
        else if (xw == 1)
        {
          this.A1 = 0;
          this.A2 = 1;
        }

        // this is the critical value of x where the used function changes
        this.xc = 1 - (this.A1/2 * Math.pow(xw,2) + xw);

        return this;
    }

    magnification(d){
        if (!arguments.length) return this.d;
        d = +d;
        if (d <= 0)
            d = 0;

        this.d = d;
        return this.rescale();
    }

    demagnificationWidth(xw){
        if (!arguments.length) return this.xw;
        xw = +xw;
        this.xw = xw;
        return this.rescale();
    }

    rangeX(_){
        if (!arguments.length) return this.range_x;
        this.range_x = _;
        return this;
    }

    range(_){
        if (!arguments.length) return this.range_x;
        return this.rangeX(_);
    }
    
    rangeY(_){
        if (!arguments.length) return this.range_y;
        this.range_y = _;
        return this;
    }
    
    focus(_){
        if (!arguments.length) return this.focus_x;
        this.focus_x = _;
        return this;
    }

    radius(_){
        if (!arguments.length) return this.r;
        this.r = +_;
        return this;
    }

    functionMode(_){
        if (!arguments.length) return this.mode;
        this.mode = _;
        return this;
    }

    fisheyeRadial(pos,inverse = false){

        if ((this.focus_x === null)|| (typeof this.focus_x == 'undefined'))
            return pos;
      
        if ((this.xw == 1 || this.d == 0)&& this.mode == 'continuous')
            return pos;

        let x = pos[0];
        let y = pos[1];

        let fx = this.focus_x[0];
        let fy = this.focus_x[1];

        let dr = Math.sqrt( Math.pow(x-fx,2) + Math.pow(y-fy,2) );

        if ((Math.abs(dr) > this.r) || (dr == 0))
        {
            //console.log('identity!');
            return pos;
        }

        let theta = Math.atan2(y-fy, x-fx);
        let cos = Math.cos(theta);
        let sin = Math.sin(theta);
        
        let dmax = this.r;
        let d = this.d;
        let xc = this.xc;
        let A1 = this.A1;
        let A2 = this.A2;

        let rescaled = dr / dmax;
        let newX = [0,0];
        let new_r;

        if (this.mode == 'continuous')
        {
            if (!inverse)
              new_r = this.fisheyeContinuous(rescaled);
            else
              new_r = this.fisheyeContinuousInverse(rescaled);
        }
        else if (this.mode == 'sarkar')
        {
            if (!inverse)
              new_r = this.fisheyeSarkar(rescaled);
            else
              new_r = this.fisheyeSarkarInverse(rescaled);
        }

        newX[0] = fx + cos * dmax * new_r;
        newX[1] = fy + sin * dmax * new_r;

        return newX;
    }

    fisheyeRadialInverse(x) {
        return this.fisheyeRadial(x, true);
    }

    fisheyeContinuous(x){
        if (x <= this.xc)            
            return (x * (this.d+1)) / (this.d*x + this.A2);
        else
            return 1 - ( -1/this.A1 + Math.sqrt(1/(this.A1*this.A1) + 2*(1-x)/this.A1) ) ;
    }
  
    fisheyeContinuousInverse(x){
        if (x <= 1 - this.xw)            
            return this.A2 * x / (this.d * (1-x) + 1);
        else
            return x - this.A1/2 * Math.pow(1-x,2) ;
    }

    fisheyeSarkar(x) {
        return (this.d+1) * x / (this.d*x + 1);
    }
  
    fisheyeSarkarInverse(x) {
        return 1 - (this.d+1) * (1-x) / (this.d*(1-x) + 1);
    }
}
