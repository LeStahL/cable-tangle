/* Cable Tangle
 * Copyright (C) 2020 Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#version 130

out vec4 gl_FragColor;

const vec2 iResolution = vec2(1920, 1080);
const float fsaa = 4.;

const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);
const float box_size = .45;

// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.
// See https://www.shadertoy.com/view/4djSRW

void hash12(in vec2 p, out float d)
{
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    d = fract((p3.x + p3.y) * p3.z);
}

void hash11(in float p, out float d)
{
    p = fract(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    d = fract(p);
}
// End of CC-A-SA 4.0 Public License

void lf11(in float t, out float n)
{
    float i = floor(t);
    t = fract(t);
    vec2 r;
    hash11(1.e2*i, r.x);
    hash11(1.e2*(i+1.), r.y);
    n = mix(r.x,r.y,smoothstep(0.,1.,t));
    n = -1.+2.*n;
}

void cartesian_to_polar(in vec2 x, out vec2 y)
{
    y = vec2(length(x), atan(x.y,x.x));
}

void polar_to_cartesian(in vec2 x, out vec2 y)
{
    y = x.x*vec2(cos(x.y),sin(x.y));
}

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = (sda.x<sdb.x)?sda:sdb;
}


float sm(in float d)
{
    return smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d);
}

void main_scene(in vec3 x, out vec2 sdf)
{
    
    const float box_size = .4;
    float r,b,ra;
    
    // Floor
	sdf = vec2(x.z+box_size,0.);    
    
    lf11(2.*x.x,b); 

	x.yz *= 1.+.4*b;

    
    float ba = 2.*x.x-b;
    mat2 R = mat2(cos(ba), sin(ba), -sin(ba), cos(ba));
    x.yz = R * x.yz;
	   
    vec2 rp,
        size = vec2(.02,pi/16.),
        y;
    cartesian_to_polar(x.yz, rp);
	vec2 drp = mod(rp, size)-.5*size,
		rpi = rp-drp;    
   	polar_to_cartesian(rpi,y);
    
    hash12(rpi, r);
    hash12(rpi+1337., ra);
    
    if(rpi.x < .3 && ra < .2)
    {
        vec2 bb;
        lf11(12.*x.x+1337.*r,bb.x); 
        lf11(12.*x.x+2337.*r,bb.y); 
        add(sdf, vec2(length((x.yz-y-.005*bb))-.02*(.125),  2.+floor(6.*r)), sdf);
    }
    
    lf11(2.*x.x-1337.,b); 
    ba = 2.*x.x-b;
    R = mat2(cos(ba), -sin(ba), sin(ba), cos(ba));
    x.yz = R*R*x.yz;
    x.x += 1337.;
    cartesian_to_polar(x.yz, rp);
	drp = mod(rp, size)-.5*size;
		rpi = rp-drp;    
   	polar_to_cartesian(rpi,y);
    
    hash12(rpi, r);
    hash12(rpi+1337., ra);
    
    if(rpi.x < .2 && ra < .2)
    {
        vec2 bb;
        lf11(12.*x.x+1337.*r,bb.x); 
        lf11(12.*x.x+2337.*r,bb.y); 
        add(sdf, vec2(length((x.yz-y-.005*bb))-.02*(.125),  2.+floor(6.*r)), sdf);
    }
    //sdf.x = abs(sdf.x)-.01;
}


float cond(in vec3 x)
{
    float b;
    lf11(2.*x.x,b); 
    x.yz *= 1.+.4*b;
    float ba = 2.*x.x-b;
    mat2 R = mat2(cos(ba), sin(ba), -sin(ba), cos(ba));
    x.yz = R * x.yz;
    return length(x.yz)-.4;
}

#define normal(o, t)void o(in vec3 x, out vec3 n, in float dx){vec2 s, na;t(x,s);t(x+dx*c.xyy, na);n.x = na.x;t(x+dx*c.yxy, na);n.y = na.x;t(x+dx*c.yyx, na);n.z = na.x;n = normalize(n-s.x);} 
normal(main_normal, main_scene)

#define march(id, sid, exit, step)void id(out vec3 x, in vec3 o, inout float d, in vec3 dir, in int N, out int i, out vec2 s){for(i = 0; i<N; ++i){x=o+d*dir;sid(x,s);if(s.x < 1.e-4) return; if(exit){i=N;}d+=step;}}
march(march_main, main_scene, false, min(s.x,mix(2.e-1,1.e-3,step(cond(x),0.))))
march(march_shadow, main_scene, cond(x)>0., min(s.x,1.e-2))

void illuminate(in vec3 x, in vec3 n, in vec3 dir, in vec3 l, inout vec3 col, in vec2 s)
{
	if(s.y == 0.) // Floor
    {
        col = vec3(.05,.08,.15);
        col = .1*col
            + .8*col*max(dot(l-x,n),0.)
            + .5*col*pow(max(dot(reflect(l-x,n),dir),0.),2.);
    }
    else if(s.y == 1.) // Structure wireframe white
    {
        col = .4*c.xxx;
        col = .6*col
            + .8*col*max(dot(l-x,n),0.)
            + 1.5*col*pow(max(dot(reflect(l-x,n),dir),0.),2.);
    }
    else if(s.y == 2.)
    {
        col = vec3(0.18,0.14,0.19);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
    else if(s.y == 3.)
    {
        col = vec3(0.29,0.27,0.47);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
    else if(s.y == 4.)
    {
        col = vec3(0.62,0.44,0.62);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
    else if(s.y == 5.)
    {
        col = vec3(0.87,0.24,0.59);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
    else if(s.y == 6.)
    {
        col = vec3(0.72,0.17,0.42);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
    else if(s.y == 7.)
    {
        col = vec3(0.59,0.07,0.29);
        col = .7*col
            + .8*col*max(dot(l-x,n),0.)
            + .3*col*pow(max(dot(reflect(l-x,n),dir),0.),1.);
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord-.5*iResolution.xy)/iResolution.y,
        s, ss,
        s0;
    vec3 col = c.yyy,
        o = 2.*c.zzx-.4*c.yyx+.05*c.zzy,
        o0 = o,
        r = c.xzy,
        t = c.yyy, 
        u = cross(normalize(t-o),-r),
        dir,
        n, 
        x,
        c1 = c.yyy,
        l,
        l2,
        dir0,
        x0,
        c2,
        n0;
    int N = 4000,
        i;
    float d = 0., d0;
    
    t += uv.x * r + uv.y * u;
    dir = normalize(t-o);

    //march_pre(x,o,d,dir,100,i,s);
	march_main(x, o, d, dir, N, i, s);
    
    l = c.yyx;
    l2 = c.xzx;
    
    if(i<N)
    {
        main_normal(x, n, 5.e-5);
        n = round(n);
    }
    else
    {
    	d = -(o.z+box_size)/dir.z;
    	x = o + d * dir;
        main_scene(x,s);
        main_normal(x, n, 5.e-5);
        n = round(n);
    }
    
    illuminate(x, n, dir, l, col, s);
    illuminate(x, n, dir, l2, c1, s);
    col = mix(col, c1, .5);
    
    if(s.y == 0.)
    {
        o0 = x;
        dir0 = reflect(dir, n);
        d0 = .01;
        march_main(x0, o0, d0, dir0, N, i, s0);
        
        if(i<N)
        {
            main_normal(x0, n0, 5.e-4);
            n = round(n);
            illuminate(x0, n0, dir0, l, c1, s0);
            illuminate(x0, n0, dir0, l2, c2, s0);
            c1 = mix(c1, c2, .5);
        	col = mix(col, c1, .4);
        }
    }

    // Soft shadow
    x0 = x;

    o = x;
    dir = normalize(l-x);
    d = 1.e-2;
    // analytical_box(o, dir, box_size*c.xxx, d);
    
    // if(d < 1.e2)
    {
        float res = 1.0;
        float ph = 1.e20;
        for(int i=0; i<N; ++i)
        // for(d=1.e-2; x.z<.5; )
        {
            x = o + d * dir;
            main_scene(x, s);
            if(s.x < 1.e-4) 
            {
                res = 0.;
                break;
            }
            if(x.z>.5) break;
            float y = s.x*s.x/(2.0*ph);
            float da = sqrt(s.x*s.x-y*y);
            res = min( res, 100.0*da/max(0.0,d-y) );
            ph = s.x;
            d += min(s.x,5.e-3);
        }
        col = mix(.3*col, col, res);
    }
    
    x = x0;
    o = x;
    dir = normalize(l2-x);
    d = 1.e-2;
    // analytical_box(o, dir, box_size*c.xxx, d);
    
    // if(d < 1.e2)
    {
        float res = 1.0;
        float ph = 1.e20;
        for(int i=0; i<N; ++i)
        // for(d=1.e-2; x.z<.5; )
        {
            x = o + d * dir;
            main_scene(x, s);
            if(s.x < 1.e-4) 
            {
                res = 0.;
                break;
            }
            if(x.z>.5) break;
            float y = s.x*s.x/(2.0*ph);
            float da = sqrt(s.x*s.x-y*y);
            res = min( res, 100.0*da/max(0.0,d-y) );
            ph = s.x;
            d += min(s.x,5.e-3);
        }
        col = mix(.3*col, col, res);
    }
    
    // Gamma
    col = col + 1.*col*col*col;
    // col *= col;
    //col = col + 1.*col*col;
    
    fragColor = vec4(clamp(col,0.,1.),1.);
}

void main()
{
    vec4 col = vec4(0.);
    float bound = sqrt(fsaa)-1.;
   	for(float i = -.5*bound; i<=.5*bound; i+=1.)
        for(float j=-.5*bound; j<=.5*bound; j+=1.)
        {
            vec4 c1;
            mainImage(c1, gl_FragCoord.xy+vec2(i,j)*1./max(bound, 1.));
     		col += c1;
        }
    col /= fsaa;
    gl_FragColor = col;
}
