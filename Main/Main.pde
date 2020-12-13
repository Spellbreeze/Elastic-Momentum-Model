// 3451-2020-P5 
// Student: Alexander Goebel
// Using base-code provided by Jarek ROSSIGNAC


import processing.pdf.*;    // to save screen shots as PDFs, does not always work: accuracy problems, stops drawing or messes up some curves !!!
import java.awt.Toolkit;
import java.awt.datatransfer.*;

import java.util.Random;

Boolean PickedFocus=false;
Boolean showFloor=true; // to show/hide the floor
  

  
void setup() 
  {
  //size(1200, 1200, P3D); // P3D means that we will do 3D graphics
  size(1000, 1000, P3D); // P3D means that we will do 3D graphics
  myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  textureMode(NORMAL);          
  noSmooth();
  
  //******************************************************** ADJUST DOWN FOR SPEED (IF NEEDED) 
  sphereDetail(24);
  
  //******************************************************** FOR SLOW MO: CHANGE THIS TO frameRate(30);  
  frameRate(30);
  
  //******************************************************** CALLS YOUR INITIALIZATION CODE FOR THIS PROJECT
  declareAllBalls();
  reinitialize();
  tc=predict();
 
  }

void draw() {
  background(255);
  hint(ENABLE_DEPTH_TEST); 

  pushMatrix();   // to ensure that we can restore the standard view before writing on the canvas
    setView();  // see pick tab
    if(showFloor) showFloor(); // draws dance floor as yellow mat
    doPick(); // sets Of and axes for 3D GUI (see pick Tab)
  
  
    //******************************************************** CALLS YOUR CODE FOR THIS PROJECT
    display(); // display balls and scaled velocity vectors
    if(animating) advanceToNextFrame(); // advances animation (processing interiediate collisions) to next frame
   
  popMatrix(); // done with 3D drawing. Restore front view for writing text on canvas

  hint(DISABLE_DEPTH_TEST); // no z-buffer test to ensure that help text is visible
  //scribeHeader(n+" balls, "+", r = "+nf(r,1,2)+ic+"ic="+ic+", jc="+jc+", tc = "+nf(tc,1,2)+", tf = "+nf(tf,1,2)+" : Event = "+Event+". Gravity = "+gravity,1);
  scribeHeader(n+" balls, "+", r="+nf(r,1,2)+", Gravity = "+gravity+". "+ 
              //"Uniform = "+uniform+", "+", tests = "+resetTestsCount+". "+
              "Event="+Event+": ic="+ic+", jc="+jc+", tc="+nf(tc,1,2)+", tf="+nf(tf,1,2)
               +", Init. Method " + ((method) ? "Spherical" : "Gaussian"),1);
  // used for demos to show red circle when mouse/key is pressed and what key (disk may be hidden by the 3D model)
  if(mousePressed) {stroke(cyan); strokeWeight(3); noFill(); ellipse(mouseX,mouseY,20,20); strokeWeight(1);}
  if(keyPressed) {stroke(red); fill(white); ellipse(mouseX+14,mouseY+20,26,26); fill(red); text(key,mouseX-5+14,mouseY+4+20); strokeWeight(1); }
  if(scribeText) {fill(black); displayHeader();} // dispalys header on canvas, including my face
  if(scribeText && !filming) displayFooter(); // shows menu at bottom, only if not filming
  if(filming && change) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif");  // save next frame to make a movie
  change=false; // to avoid capturing frames when nothing happens (change is set uppn action)
  change=true;
  }
