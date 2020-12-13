//******************************************************** PROJECT CONTROLS
boolean snapping; // make a picture of each display frame
boolean slowMo; // toggle slow motion (1 fpg or 30 fps)

boolean phaseA = false;
boolean phaseB = false;
boolean phaseC = false;

//******************************************************** CONTAINER
float rc = 600; // radius of container
PNT C = P(0,0,0); // container

//******************************************************** BALLS
int maxn = 4096; //TODO: reduce for Phases 1&3
int n=1;                 // current number of balls ('<','>')
PNT[] B = new PNT[maxn]; // centers of ball
float r = 30;           // radius of all balls

//******************************************************** VELOCITIES
VCT[] V = new VCT[maxn]; // velocities of balls
boolean showVelocities=true;
float m = 2; // initial speed of balls (magnitude of V[])
float s = 10;  // scaling of arrows when displaying vectors V[]

//******************************************************** FLAT VS 3D
float z=0; // height: when z==0  ball centers stay on to the floor
boolean flat=true; // true when ball centers stay on the floor

//******************************************************** GRAVITY
boolean gravity=false; // toggle to add gravity
VCT Gravity=V(0,0.4,0); // gravity vector  for 2D

//******************************************************** ANIMATION CONTROL
boolean animating=false;            // automatic animation mode
boolean advancingToNextFrame=false;  // boolean set by 'a'
float dtf = 1;                   // interframe time-lapse (1/30 sec) is one unit of time
float tc = 0;                      // remainig time to next collision
float tf = dtf;                      // remainig time to next frame
String Event ="";
int ic=-1, jc=-1;                 // IDs of colliders (-1 means not collider)

//******************************************************** COLORING BALLS TO VISUALIZE EVENTS
boolean changingColors=true; // to show which balls will collide
boolean[] X = new boolean[maxn];  // mark balls that interfere with other balls or stick out of the container (for validation) 


//******************************************************** DECLARE BALLS, VELOCITIES, ATTRIBUTES
void declareAllBalls() 
  {
  for(int i=0; i<maxn; i++) B[i]=P();
  for(int i=0; i<maxn; i++) V[i]=V();
  for(int i=0; i<maxn; i++) X[i]=false; // overlapping of escaping balls
  }  


//******************************************************** INITIALIZE BALLS AND VELOCITIES
void reinitialize() // reset B[0]...B[n] to uniformely random sampling of non-overlapping balls and sets random velocities V[i] 
  {
  for(int i=0; i<maxn; i++) X[i]=false;
  //... FIx me to ensure that balls are in the container and disjoint
  int out_count = 0, interfere_count = 0;
  for(int i=0; i<n; i++) {
    while(!X[i]) {
      B[i]= RandomPoint();
      X[i] = true;
      for(int j=0; j < i; j++) {
        if(interfere(i,j)) {
          X[i] = false;
          interfere_count++;
          break;
        }
      }
      
      if (out_count + interfere_count > 80000) {
        System.out.println("Reverted to n/2 - Error Count: " + out_count + interfere_count);
        n /= 2;
        i = -1;
        out_count = 0; interfere_count = 0;
        break;
      }
    }
  }
  System.out.println("Error Count for " + n + " balls- Interfering: " + interfere_count + "\t Out of Region: " + out_count);
  check(); // to test whether balls overlap
  for(int i=0; i<n; i++) {
    if(flat) {float w = random(TWO_PI); V[i]=V(m*cos(w),m*sin(w),0);}
    else V[i]=V(m,RandomDirection());
  }
 }

boolean case2_3 = false;
void fixedReinitialize() { //For testing purposes
  if (case2_3) {
      if (method) {
        n = 4;
        B[0] = P(0,3*r,0);
        B[1] = P(0,r,0);
        B[2] = P(0,-r,0);
        B[3] = P(0,-6*r,0);
        V[0] = V();
        V[1] = V();
        V[2] = V();
        V[3] = V(0,18,0);
      } else {
        n = 4;
        B[0] = P(0,6*r,0);
        B[1] = P(0,r,0);
        B[2] = P(0,-r,0);
        B[3] = P(0,-6*r,0);
        V[0] = V(0,-18,0);
        V[1] = V();
        V[2] = V();
        V[3] = V(0,18,0);
      }
  } else if (method) {
      n = 3;
      B[0] = P(0,-3*r,0);
      B[1] = P(3*r*sqrt(3)/2,3*r/2,0);
      B[2] = P(-3*r*sqrt(3)/2,3*r/2,0);
      V[0] = V(0,3,0);
      V[1] = V(-3.0*sqrt(3)/2, -3.0/2,0);
      V[2] = V(3.0*sqrt(3)/2, -3.0/2,0);
      V[3] = V();
      V[4] = V();
      V[5] = V();
    } else {
      n = 4;
      B[0] = P(-3*r,-4*r,0);
      B[1] = P(0,-r,0);
      B[2] = P(0,r,0);
      B[3] = P(-3*r,4*r,0);
      //B[5] = P(3*r,2*r,0);
      V[0] = V(1,1,0);
      V[1] = V();
      V[2] = V();
      V[3] = V(1,-1,0);
      V[4] = V();
      V[5] = V();
    }
}

//******************************************************** INTERFERENCE WITH CONTAINER AND COLLISION TESTS
boolean interfere(int i, int j) // balls i & j interfere 
  {
  return V(B[i],B[j]).norm() < 2*r;
  }  

boolean sticksOut(int i) // balls interferes with container 
  {
  return V(C, B[i]).norm() > rc-r;
  }  

//******************************************************** UNIFORMLY RANDOM DIRECTIONS IN 3D
VCT RandomDirection() 
  {
  //... Fix me to produce uniformly distributed random directions
  return U(gauss(),gauss(),z*gauss()); // z==0 when flat (planar motion)
  }
  
//gives a uniformly random point
boolean method = true;
PNT RandomPoint() {
  if (method) { // this method seems to cluster around the axes
    VCT point = V(1,0,0);
    float angle_y = rand()*360-180;
    float angle_z = rand()*360-180;
    return P(C,V((float) Math.pow(rand(),1.0/3) * (rc-r), U((point.rotate(angle_y, V(0,0,1)).rotate(z*angle_z,V(0,1,0))))));
  } else { //https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability/87238#87238
    VCT point = U(V(gauss(), gauss(), z*gauss()));
    point = V((rc-r) * (float) Math.pow(rand(), 1.0/3), point);
    return P(C,point);
  }
}

float rand() {
  return (float) Math.random();
}
Random gaussPicker = new Random();
float gauss() {
  return (float) gaussPicker.nextGaussian(); 
}
  
void scaleVelocties(float pm) {for(int i=0; i<n; i++) V[i].mul(pm);}

void check() // for testing only: sets X[i]=true; X[j]=true; if balls i and j interfere
  {
  //...
  }
  

//******************************************************** DISPLAY ALL COLORED BALLS AND VELOCITIES
void display() // shows balls and velocities and colors colliders
  {
  fill(blue); if(showVelocities) for(int i=0; i<n; i++) arrow(B[i],s,V[i],10);  
  fill(cyan);
  for(int i=0; i<n; i++) 
    {
    if(changingColors)
      if(i==ic || i==jc) {if(jc==-1) fill(magenta); else fill(red);} 
      else {if(X[i]) fill(dgreen); else fill(cyan);} 
    show(B[i],r);
    } 
  fill(orange,100); for(int i=0; i<n; i++) if(X[i]) show(B[i],r+5); 
  fill(color(150,200,100,100)); 
  pushMatrix(); if(flat) scale(1,1,0.01); // flat mode squashes the container
  show(C,rc);
  popMatrix();
  }  

 
//******************************************************** ADVANCES ANMATION TO NECT DISPLAY-FRAME 
void advanceToNextFrame() // may involve any number of collisions or bounces along the way
  {
  for (int i=0; i<n; i++) {
    if (gravity) {
       V[i] = A(V[i],Gravity); 
    }
    show(B[i], r); 
  }
  if (!phaseA) {   
    tc=predict();
    boolean more_collisions = false;
    do {
      if (tc < tf) {
         tf -= tc;
         jump(tc);
         clash();
         tc=predict();
         more_collisions = tc <= tf;
      } else {
        jump(tf);
        tf = dtf;
      }
    } while (more_collisions);
  } else {
    // PHASE A Code
    tc=predict();
    float tt = 0;
    float td = tf/20; // time delta
    float dtd = td; // default time delta
    float toac = 0; //time of last collision
    //boolean first_collision = true;
    boolean more_collisions;
    boolean coll_occur;
    while (tt<=tf) {
      more_collisions = false;
      coll_occur = false;
      do {
        if (tc < td) {
          td -= tc;
          jump(tc);
          clash();
          tc=predict();
          coll_occur = true;
          more_collisions = tc <= td;
        } else {
          jump(td);
          tc -= td;
          td = dtd;
        }
      } while (more_collisions);
      if (coll_occur) { 
        fill(cyan);
      } else {
        fill(250,250,0,80); //yellow (not bouncing);
      }
      for (int i=0; i<n; i++) {
        show(B[i], r); 
      }
      tt += td;
    }
    fill(red);
    for(int i=0; i<n; i++) {
      show(B[i], r+1); 
    }
  }
  
  }

void jump(float t) {
  for(int i=0; i<n; i++) B[i]=P(B[i],t,V[i]);
} // advances all balls by time t

void advanceToPreCollision() {
  while (tc > tf) {
    jump(tf); 
    tc=predict();
  }
}
//******************************************************** PREDICT NEXT COLLISION/BOUNCE
float predict() // compute remainingTimeToFirstClash and colliders' IDs (ic,jc)
  {
    float ltc=100000, edgeBounceTime, collisionBounceTime;
    ic=-1; jc=-1;
    for(int i=0; i<n; i++) {
      edgeBounceTime = timeToBounce(i);
      if (edgeBounceTime < ltc && edgeBounceTime >= 0) {
        ltc = edgeBounceTime;
        ic = i;
        jc = -1;
      }
      for(int j=0; j<i; j++) {
        collisionBounceTime = timeToCollision(i,j);
        if (collisionBounceTime < ltc && collisionBounceTime >= 0) {
          ltc = collisionBounceTime;
          ic = i;
          jc = j;
        }
      }
    }
    return ltc;
  }  
  
float timeToCollision(int i, int j) // time to collisions between balls i and j, or -1 if no collision
  {   
    VCT w_minus_v = M(V[j], V[i]);
    VCT BA = V(B[i], B[j]);
    
    // (t^2)(W-V)dot(W-v) + (t)(2*((W-V)dot(BA))) + ((BA)dot(BA) - 4r^2))
    float pyth_A = dot(w_minus_v,w_minus_v); // velocity dotted with itself
    float pyth_B = 2*dot(w_minus_v,BA);
    float pyth_C = dot(BA,BA)-4*r*r;
    float plus_t = (float) (-pyth_B + Math.pow(pyth_B*pyth_B - 4*pyth_A*pyth_C,0.5))/(2*pyth_A);
    float minus_t = (float) (-pyth_B - Math.pow(pyth_B*pyth_B - 4*pyth_A*pyth_C,0.5))/(2*pyth_A);
    
    if (minus_t < 0) {// to mitigate the floating-point error
      minus_t = (minus_t > -0.000001) ? 0 : minus_t;
    }
    if (plus_t < 0) {// to mitigate the floating-point error
      plus_t = (plus_t > -0.00001) ? 0 : plus_t;
    }
    
    float retval = -1;
    if (plus_t >= 0) {
      if (minus_t >= 0) {
        retval = Math.min(plus_t,minus_t);
      } else {
        retval = plus_t;
      }
    } else if (minus_t >= 0) {
      retval = minus_t;
    }
    
    if (retval != -1 && approaching(i,j,retval)) {
       return retval; 
    } else {
      return -1;
    }//*/
  } 

float timeToBounce(int i) // time to bounce with collider's border for ball i , or -1 if no collision
  {
    VCT CB = V(C, B[i]);
    //(t^2)((V)dot(V)) + (t)(2*(CB)dot(V)) + ((CB)dot(CB) - (rc - r)^2)
    float pyth_A = dot(V[i],V[i]);
    float pyth_B = 2 * dot(CB, V[i]);
    float pyth_C = dot(CB,CB) - sq(rc-r);
    
    float plus_t = (float) (-pyth_B + Math.pow(pyth_B*pyth_B - 4*pyth_A*pyth_C,0.5))/(2*pyth_A);
    float minus_t = (float) (-pyth_B - Math.pow(pyth_B*pyth_B - 4*pyth_A*pyth_C,0.5))/(2*pyth_A);
    
    if (plus_t < 0) {// to mitigate the floating-point error
      plus_t = (plus_t > -0.00001) ? 0 : plus_t;
    }
    
    
    float retval = -1;
    
    return (plus_t >= 0) ? plus_t : -1; //Only plus_t because always inside of sphere
  }    

//******************************************************** TEST WHETHER BALL(S) MOVING TOWARDS COLLISION
boolean approaching(int i, int j, float t) // balls i & j will be aproaching each other  at time t
  {
    return d(P(B[i],t+0.001,V[i]),P(B[j],t+0.001,V[j]))
              < d(P(B[i],t,V[i]),P(B[j],t,V[j]));
  }  

boolean exiting(int i, float t) // balls i will overlap with the container border from inside at time t
  {
    return V(C, P(B[i],t,V[i])).norm() >= (rc-r) && V(C, P(B[i],t,V[i])).norm() < V(C, P(B[i],t+0.001,V[i])).norm();
  }  



//******************************************************** COMPUTE NEW VELOCITIES AFTER COLLISION OR BOUNCE
void clash() // computes new velocities after elasic collision with container or other ball
  {
  if(jc==-1) {computeNewVelocitiesAfterBounceOffContainer(ic);}
  else {
    computeNewVelocitiesAfterBounceBetweenColliders(ic,jc);
    }
    ic=-1;
    jc=-1;
  }    

void computeNewVelocitiesAfterBounceOffContainer(int i) 
  {
  VCT new_veloc = V[i];
  PNT bounce_point = P(C, A(V(C,B[i]),V(r, U(V(C,B[i])))));
  VCT normal;
  normal = V(B[i], bounce_point).normalize();
  normal = V(dot(new_veloc,normal),normal);
  new_veloc = M(new_veloc, normal);
  
  if (gravity) {
     normal = V(0.99, normal);
    }
  V[i] = M(new_veloc,normal);
  if (gravity) {
       if((projection(V[ic],Gravity)).norm() > rc*0.5) {
        V[ic] = V(); 
       }
    }  
  }

void computeNewVelocitiesAfterBounceBetweenColliders(int i, int j)
  {
  VCT new_veloc_I = V[i], new_veloc_J = V[j];
  VCT normal_on_I, normal_on_J;
  normal_on_I = projection(V[i], V(B[i], B[j]).normalize()); //projection of Vi on normal of IJ
  normal_on_J = projection(V[j], V(B[j], B[i]).normalize()); //projection of Vj on normal of JI
  
  new_veloc_I = M(new_veloc_I, normal_on_I);
  new_veloc_J = M(new_veloc_J, normal_on_J);
  
  if (gravity) {
     normal_on_I = V(0.99, normal_on_I);
     normal_on_J = V(0.99, normal_on_J);
  }
  
  V[i] = A(new_veloc_I, normal_on_J);
  V[j] = A(new_veloc_J, normal_on_I);
  
  }
  
VCT projection(VCT V, VCT U) { //((V*U)/(U*U))*U
  return V((dot(V,U)/dot(U,U)), U);
}


//Gravity Quartic Equation https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/auckly29.pdf
/*float[] gravQuartic(float a, float b, float c, float d) {
   float p = b-(3*a*a/8);
   float q = c - a*b*0.5 + a*a*a/8;
   float r = d - a*c*0.25 + a*a*b/16 - 3*a*a*a*a/256;
   
}
float lamba(float p, float q, float r) {
  Cubic equ_cubic = new Cubic();
}*/
