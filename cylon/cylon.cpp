/**
 * Import the necessary libraries
 */
#include "fstream"
#include <sstream>
#include <string>
#include <iostream>

#include "math.h"
#include "vector"
#include "list"

#if defined(_MSC_VER)
#include <gl/glut.h>
#else
#include <GLUT/glut.h>
#endif

using namespace std;

/**
 * Define constants
 */
#define array vector // to avoid name collions
#define forXYZ(i) for(int i=0; i<3; i++)
const float PRECISION = 0.0001;
const int X=0;
const int Y=1;
const int Z=2;
const int W=3; // for quaternions only
const float gravity = 9.8; // meters/second^2

/**
 * Define class Vector
 */
class Vector {
public:
  float v[3];
  Vector(float x=0, float y=0, float z=0) {
    v[X]=x; v[Y]=y; v[Z]=z;
  }
  float operator()(int i) const { return v[i]; }
  float &operator()(int i) { return v[i]; }
};

/**
 * Define operations between vectors
 */
Vector operator*(float c, const Vector &v) {
  return Vector(c*v(X),c*v(Y),c*v(Z));
}
Vector operator*(const Vector &v, float c) {
  return c*v;
}
Vector operator/(const Vector &v, float c) {
  return (1.0/c)*v;
}
Vector operator+(const Vector &v, const Vector &w) {
  return Vector(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z));
}
Vector operator-(const Vector &v, const Vector &w) {
  return Vector(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z));
}
float operator*(const Vector &v, const Vector &w) {
  return v(X)*w(X) + v(Y)*w(Y) + v(Z)*w(Z);
}
Vector cross(const Vector &v, const Vector &w) {
  return Vector(v(Y)*w(Z)-v(Z)*w(Y),
		 v(Z)*w(X)-v(X)*w(Z),
		 v(X)*w(Y)-v(Y)*w(Z));
}

/**
 * Class Rotation
 */
class Matrix {
public:
  float m[3][3];
  Matrix() { forXYZ(i) forXYZ(j) m[i][j]=0; }
  const float operator()(int i, int j) const { return m[i][j]; }
  float &operator()(int i, int j) { return m[i][j]; }
};

class Rotation : public Matrix {
public:
  Rotation() { m[X][X]=m[Y][Y]=m[Z][Z]=1; }
  Rotation(const Vector& v) {
    float theta = sqrt(v*v);
    if(theta<PRECISION) {
      m[X][X]=m[Y][Y]=m[Z][Z]=1;
    } else {
      float s = sin(theta), c=cos(theta);
      float t = 1-c;
      float x = v(X)/theta, y = v(Y)/theta, z = v(Z)/theta;
      m[X][X]=t*x*x+c;   m[X][Y]=t*x*y+s*z; m[X][Z]=t*x*z-s*y;
      m[Y][X]=t*x*y-s*z; m[Y][Y]=t*y*y+c;   m[Y][Z]=t*y*z+s*x;
      m[Z][X]=t*x*z+s*y; m[Z][Y]=t*y*z-s*x; m[Z][Z]=t*z*z+x;
    }
  }
};

/**
 * Class Inertia Tensor
 */
class InertiaTensor : public Matrix {};
InertiaTensor operator+(const InertiaTensor &a, const InertiaTensor &b) {
  InertiaTensor c=a;
  forXYZ(i) forXYZ(j) c(i,j)+=b(i,j);
  return c;
}

/**
 * Operations between matrices
 */
Vector operator*(const Matrix &R, const Vector &v) {
  return Vector(R(X,X)*v(X)+R(X,Y)*v(Y)+R(X,Z)*v(Z),
		 R(Y,X)*v(X)+R(Y,Y)*v(Y)+R(Y,Z)*v(Z),
		 R(Z,X)*v(X)+R(Z,Y)*v(Y)+R(Z,Z)*v(Z));
}
Matrix operator*(const Matrix &R, const Matrix &S) {
  Matrix T;
  forXYZ(i) forXYZ(j) forXYZ(k) T(i,j)+=R(i,k)*S(k,j);
  return T;
}
float det(const Matrix &R) {
  return R(X,X)*(R(Z,Z)*R(Y,Y)-R(Z,Y)*R(Y,Z))
    -R(Y,X)*(R(Z,Z)*R(X,Y)-R(Z,Y)*R(X,Z))
    +R(Z,X)*(R(Y,Z)*R(X,Y)-R(Y,Y)*R(X,Z));
}
Matrix operator/(float c, const Matrix &R) {
  Matrix T;
  float d = c/det(R);
  T(X,X)=(R(Z,Z)*R(Y,Y)-R(Z,Y)*R(Y,Z))*d;
  T(X,Y)=(R(Z,Y)*R(X,Z)-R(Z,Z)*R(X,Y))*d;
  T(X,Z)=(R(Y,Z)*R(X,Y)-R(Y,Y)*R(X,Z))*d;
  T(Y,X)=(R(Z,X)*R(Y,Z)-R(Z,Z)*R(Y,X))*d;
  T(Y,Y)=(R(Z,Z)*R(X,X)-R(Z,X)*R(X,Z))*d;
  T(Y,Z)=(R(Y,X)*R(X,Z)-R(Y,Z)*R(X,X))*d;
  T(Z,X)=(R(Z,Y)*R(Y,X)-R(Z,X)*R(Y,Y))*d;
  T(Z,Y)=(R(Z,Y)*R(X,Y)-R(Z,Y)*R(X,X))*d;
  T(Z,Z)=(R(Y,Y)*R(X,X)-R(Y,X)*R(X,Y))*d;
  return T;
}

/**
 * Class Body (describes a rigid body object)
 */
class Body {
public:
  // object shape
  array<Vector> r;  // vertices in local coordinates
  array<array<int> > faces;
  // properties of the body /////////////////////////
  bool locked;       // if set true, don't integrate
  float m;           // mass
  InertiaTensor I;   // moments of inertia (3x3 matrix)
  // state of the body //////////////////////////////
  Vector p;          // position of the center of mass
  Vector K;          // momentum
  Matrix R;          // orientation
  Vector L;          // angular momentum
  // auxiliary variables ///////////////////////////
  float inv_m;       // 1/m
  Matrix inv_I;      // 1/I
  Vector F;
  Vector tau;
  Vector v;          // velocity
  Vector omega;      // angular velocity
  array<Vector> Rr; // rotated r's.
  array<Vector> vertices; // rotated and shifted r's
  // forces and constraints
  // ...
  Body(float m=1.0, bool locked=false) {
    this->locked = locked;
    this->m = 1.0;
    I(X,X)=I(Y,Y)=I(Z,Z)=m;
    inv_m = 1.0/m;
    inv_I = 1.0/I;
    R = Rotation();
  };
  void clear() { F(X)=F(Y)=F(Z)=tau(X)=tau(Y)=tau(Z)=0; }
  void update_vertices();
  void integrator(float dt);
  void loadObj(const string & file);
};

/**
 * rotate and shift all vertices from local to universe
 */
void Body::update_vertices() {
  Rr.resize(r.size());
  vertices.resize(r.size());
  for(int i=0; i<r.size(); i++) {
     Rr[i] = R*r[i];
     vertices[i]=Rr[i]+p;
  }
}

/**
 * Euler integrator
 */
void Body::integrator(float dt) {
  v     = inv_m*K;
  omega = inv_I*L;
  p     = p + v*dt;                // shift
  K     = K + F*dt;                // push
  R     = Rotation(omega*dt)*R;    // rotate
  L     = L + tau*dt;              // spin
  update_vertices();
}

/**
 * Interface for all forces.
 * The constructor can be specific of the force
 * the apply methods adds the contribution to F and tau
 */
class Force {
public:
  virtual void apply()=0;
};

/**
 * Gravity
 */
class GravityForce : public Force {
public:  
  Body *body;
  float g;
  GravityForce(Body *body, float g = gravity) {
    this->body=body; this->g=g;    
  }
  void apply() {
    body->F(Y) -= (body->m)*g;
  }
};

class SpringForce : public Force {
public:
  Body *bodyA;
  Body *bodyB;
  int iA, iB;
  float kappa, L;
  SpringForce(Body *bodyA, int iA, Body *bodyB, int iB,
	      float kappa, float L) {
    this->bodyA=bodyA; this->iA=iA;
    this->bodyB=bodyB; this->iB=iB;
  }
  void apply() {
    Vector d = bodyA->vertices[iA]-bodyB->vertices[iB];
    float n = sqrt(d*d);
    Vector F = kappa*(L-n)*d/n;
    bodyA->F = bodyA->F+F;
    bodyB->F = bodyA->F-F;
    bodyA->tau = bodyA->tau + cross(bodyA->Rr[iA],F);
    bodyB->tau = bodyB->tau - cross(bodyB->Rr[iB],F);
  }
};

class AnchoredSpringForce : public Force {
public:
  Body *body;
  int i;
  Vector pin;
  float kappa, L;
  AnchoredSpringForce(Body *body, int i, Vector pin,
		      float kappa, float L) {
    this->body=body; this->i=i;
    this->pin = pin;
  }
  void apply() {
    Vector d = body->vertices[i]-pin;
    float n = sqrt(d*d);
    Vector F = kappa*(L-n)*d/n;
    body->F = body->F+F;
    body->tau = body->tau + cross(body->Rr[i],F);
  }
};


/**
 * Friction (ignores the shape of the body, assumes a sphere)
 */
class FrictionForce: public Force {
public:
  Body *body;
  float gamma;
    FrictionForce(Body *body, float gamma) {
      this->body=body; this->gamma=gamma;
    }
    void apply() {
      body->F = body->F-gamma*body->v; // ignores shape
    }
};

/**
 * class that stores all bodies, forces and constraints
 */ 
class Universe {
public:
  list<Body*> bodies;
  list<Force*> forces;
  list<Body*>::iterator body;
  list<Force*>::iterator force;
  // evolve universe
  void evolve(float dt) {
    // clear forces and troques
    for(body=bodies.begin(); body!=bodies.end(); body++)
      (*body)->clear();
    // compute forces and torques
    for(force=forces.begin(); force!=forces.end(); force++)
      (*force)->apply(); // adds to F and tau
    // integrate
    for(body=bodies.begin(); body!=bodies.end(); body++) 
      if(!(*body)->locked)
	(*body)->integrator(dt);
    // handle collisions
    for(body=bodies.begin(); body!=bodies.end(); body++)             
      if(!(*body)->locked) {
	if((*body)->p(Y)<0) {
	  (*body)->p(Y)=0;
	  (*body)->K(Y)=-0.9*(*body)->K(Y);
	  (*body)->L = 0.9*(*body)->L;
	}
      }
  }
};

/**
 * Auxiliary functions translate moments of Inertia
 */
 InertiaTensor dI(float m, const Vector &r) {
  InertiaTensor I;
  float r2 = r*r;
  forXYZ(j) forXYZ(k) I(j,k) = m*((j==k)?r2:0-r(j)*r(k));
  return I;
}

/**
 * Auxliary function to include to merge two bodies
 * needs some more work....
 */
Body operator+(const Body &a, const Body &b) {
  Body c;
  c.m = (a.m+b.m);
  c.p = (a.m*a.p + b.m+b.p)/c.m;
  c.K = a.K+b.K;
  c.L = a.L+b.L;
  Vector da = a.p-c.p;
  Vector db = b.p-c.p;
  c.I = a.I+dI(a.m,da)+b.I+dI(b.m,db);
  int n = a.r.size();
  // copy all r
  for(int i=0; i<n; i++)
    c.r.push_back(a.r[i]+da);
  for(int i=0; i<b.r.size(); i++)
    c.r.push_back(b.r[i]+db);
  // copy all faces and re-label r  
  int m = a.faces.size();
  c.faces.resize(a.faces.size()+b.faces.size());
  for(int j=0; j<a.faces.size(); j++)
    c.faces[j]=a.faces[j];
  for(int j=0; j<b.faces.size(); j++)
    for(int k=0; k<b.faces[j].size(); k++)
      c.faces[j+m].push_back(b.faces[j][k]+n);  
  c.update_vertices();
  return c;
}

void Body::loadObj(const string & file) {
  ifstream input;
  input.open(file.c_str());
  string line;
  if(input.is_open()) {
    while(input.good()) {
      std::getline(input, line);
      if(line.length()>0) {
	string initialVal;
	istringstream instream;
	instream.str(line);
	instream >> initialVal;
	if(initialVal=="v") {
	  float x,y,z;
	  instream >> x >> y >> z;
	  r.push_back(Vector(x,y,z));
	} else if (initialVal=="f") {
	  int v1, v2, v3;
	  instream >> v1 >> v2 >> v3;
	  array<int> triangle;
	  triangle.push_back(v1-1);
	  triangle.push_back(v2-1);
	  triangle.push_back(v3-1);
	  faces.push_back(triangle);
	}
      }
    }
    update_vertices();
  }
}


/**
 * GLUT code below
 * Creates a window in which to display the scene.
 */
void createWindow(const char* title) {
  int width = 640;
  int height = 480;
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(width,height);
  glutInitWindowPosition(0,0);
  glutCreateWindow(title);
  
  glClearColor(0.9f, 0.95f, 1.0f, 1.0f);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (double)width/(double)height, 1.0, 500.0);
  glMatrixMode(GL_MODELVIEW);
}

Universe myuniverse;
int frame=0;

/**
 * Called each frame to update the 3D scene. Delegates to
 * the application.
 */
void update() {
  // evolve world
  float timeStep = 0.016f;	// 60fps fixed rate.
  myuniverse.evolve(timeStep);
  glutPostRedisplay();
  frame+=1;
}

/**
 * Called each frame to display the 3D scene. Delegates to
 * the application.
 */
void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(0.0, 3.5, 8.0,  0.0, 3.5, 0.0,  0.0, 1.0, 0.0);
  glColor3f(0,0,0);
  
  // draw
  for(list<Body*>::iterator body=myuniverse.bodies.begin();
      body!=myuniverse.bodies.end(); body++) {
    glPushMatrix();
    Vector &pos = (*body)->p;
    glTranslatef(pos.v[X], pos.v[Y], pos.v[Z]);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    for(int i=0; i<(*body)->faces.size(); i++) {      
      glBegin(GL_POLYGON);
      for(int j=0; j<(*body)->faces[i].size(); j++) {
	int k = (*body)->faces[i][j]; 
	glVertex3fv((*body)->vertices[k].v);
      }
      int k = (*body)->faces[i][0]; 
      glVertex3fv((*body)->vertices[k].v);
      glEnd();
    }
    glPopMatrix();
  }
  // update the displayed content
  glFlush();
  glutSwapBuffers();
}

/**
 * Called when the display window changes size.
 */
void reshape(int width, int height) {
  glViewport(0, 0, width, height);
}

/**
 * Called when a mouse button is pressed. Delegates to the
 * application.
 */
void mouse(int button, int state, int x, int y) { }

/**
 * Called when a key is pressed.
 */
void keyboard(unsigned char key, int x, int y) {
  // Note we omit passing on the x and y: they are rarely needed.
  //Process keyboard
}

/**
 * Called when the mouse is dragged.
 */
void motion(int x, int y) { }

/**
 * The main entry point. We pass arguments onto GLUT.
 */
int main(int argc, char** argv) {
  // Set up GLUT and the timers
  glutInit(&argc, argv);
  
  // Create the application and its window
  createWindow("GPNS");
  Body b = Body();
  b.loadObj("assets/sphere.obj");
  b.F = Vector(0,-0.01,0);
  b.p = Vector(0,2,0);
  b.L = Vector(0.1,0,0);
  myuniverse.bodies.push_back(&b);
  myuniverse.forces.push_back(new GravityForce(&b,0.1));

  // Set up the appropriate handler functions
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutDisplayFunc(display);
  glutIdleFunc(update);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  
  // Run the application
  glutMainLoop();
  
  // Clean up the application
  return 0;
}
