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
#define forXYZ(i) for(int i=0; i<3; i++)
const int X=0;
const int Y=1;
const int Z=2;
const int W=3; // for quaternions only
const float g = 9.8; // meters/second^2

/**
 * Define class Vector3
 */
class Vector3 {
public:
  float v[3];
  Vector3(float x=0, float y=0, float z=0) {
    v[X]=x; v[Y]=y; v[Z]=z;
  }
  float operator()(int i) const { return v[i]; }
  float &operator()(int i) { return v[i]; }
};

/**
 * Define operations between vectors
 */
Vector3 operator*(float c, const Vector3 &v) {
  return Vector3(c*v(X),c*v(Y),c*v(Z));
}
Vector3 operator*(const Vector3 &v, float c) {
  return c*v;
}
Vector3 operator/(const Vector3 &v, float c) {
  return (1.0/c)*v;
}
Vector3 operator+(const Vector3 &v, const Vector3 &w) {
  return Vector3(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z));
}
Vector3 operator-(const Vector3 &v, const Vector3 &w) {
  return Vector3(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z));
}
float operator*(const Vector3 &v, const Vector3 &w) {
  return v(X)*w(X) + v(Y)*w(Y) + v(Z)*w(Z);
}
Vector3 cross(const Vector3 &v, const Vector3 &w) {
  return Vector3(v(Y)*w(Z)-v(Z)*w(Y),
		 v(Z)*w(X)-v(X)*w(Z),
		 v(X)*w(Y)-v(Y)*w(Z));
}

/**
 * Define class Quaternion
 */
class Quaternion {
public:
  float q[4];
  Quaternion(float x=0, float y=0, float z=0, float w=1) {
    q[X]=x; q[Y]=y; q[Z]=z; q[W]=w;
  }
  Quaternion(float theta, const Vector3 &v) {
    float n = sqrt(v(X)*v(X)+v(Y)*v(Y)+v(Z)*v(Z));
    float s = sin(theta/2);
    q[X]=v(X)/n*s; q[Y]=v(Y)/n*s; q[Z]=v(Z)/n;
    q[W]=cos(theta/2);
  }  
  float operator()(int i) const { return q[i]; }
  float &operator()(int i) { return q[i]; }
  void normalize() {
    float n = sqrt(q[X]*q[X]+q[Y]*q[Y]+q[Z]*q[Z]+q[W]*q[W]);
    q[X]/=n; q[Y]/=n; q[Z]/=n; q[W]/=n;
  }
};

/**
 * Opartions among quaternions
 */
Quaternion operator*(float c, const Quaternion &v) {
  return Quaternion(c*v(X),c*v(Y),c*v(Z),c*v(W));
}
Quaternion operator*(const Quaternion &v, float c) {
  return c*v;
}
Quaternion operator+(const Quaternion &v,const Quaternion &w) {
  return Quaternion(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z),v(W)+w(W));
}
Quaternion operator-(const Quaternion &v,const Quaternion &w) {
  return Quaternion(v(X)+w(X),v(Y)+w(Y),v(Z)+w(Z),v(W)-w(W));
}
Quaternion operator*(const Vector3 &v, const Quaternion &w) {
  return Quaternion(-v(X)*w(X) - v(Y)*w(Y) - v(Z)*w(Z),
		    v(X)*w(W) + v(Y)*w(Z) - v(Z)*w(Y),
		    v(Y)*w(W) + v(Z)*w(X) - v(X)*w(Z),
		    v(Z)*w(W) + v(X)*w(Y) - v(Y)*w(X));
}

/**
 * Class Matrix33
 */
class Matrix33 {
public:
  float m[3][3];
  const float operator()(int i, int j) const { return m[i][j]; }
  float &operator()(int i, int j) { return m[i][j]; }
};

/**
 * Class Rotation
 */
class Rotation : public Matrix33 {
public:
  Rotation() { m[X][X]=m[Y][Y]=m[Z][Z]=1; }
  Rotation(const Vector3& v) {
    float theta = sqrt(v*v);
    float s = sin(theta);
    float c=cos(theta);
    float t = 1-c;
    float x = v(X)/theta;
    float y = v(Y)/theta;
    float z = v(Z)/theta;
    m[X][X]=t*x*x+c;
    m[X][Y]=t*x*y+s*z;
    m[X][Z]=t*x*z-s*y;
    m[Y][X]=t*x*y-s*z;
    m[Y][Y]=t*y*y+c;
    m[Y][Z]=t*y*z+s*x;
    m[Z][X]=t*x*z+s*y;
    m[Z][Y]=t*y*z-s*x;
    m[Z][Z]=t*z*z+x;
  }
  void set(const Quaternion &theta) {
    m[X][X]=2.0*(theta(X)*theta(X)+theta(W)*theta(W))-1;
    m[X][Y]=2.0*(theta(X)*theta(X)-theta(Z)*theta(W));
    m[X][Z]=2.0*(theta(X)*theta(Z)+theta(Y)*theta(W));
    m[Y][X]=2.0*(theta(X)*theta(Y)+theta(Z)*theta(W));
    m[Y][Y]=2.0*(theta(Y)*theta(Y)+theta(W)*theta(W))-1;
    m[Y][Z]=2.0*(theta(Y)*theta(Z)-theta(X)*theta(W));
    m[Z][X]=2.0*(theta(X)*theta(Z)-theta(Y)*theta(W));
    m[Z][Y]=2.0*(theta(Y)*theta(Z)+theta(X)*theta(W));
    m[Z][Z]=2.0*(theta(Z)*theta(Z)+theta(W)*theta(W))-1;
  }
};

/**
 * Class Inertia Tensor
 */
class InertiaTensor : public Matrix33 {};
InertiaTensor operator+(const InertiaTensor &a, const InertiaTensor &b) {
  InertiaTensor c=a;
  forXYZ(i) forXYZ(j) c(i,j)+=b(i,j);
  return c;
}

/**
 * Operations between matrices
 */
Vector3 operator*(const Matrix33 &R, const Vector3 &v) {
  return Vector3(R(X,X)*v(X)+R(X,Y)*v(Y)+R(X,Z)*v(Z),
		 R(Y,X)*v(X)+R(Y,Y)*v(Y)+R(Y,Z)*v(Z),
		 R(Z,X)*v(X)+R(Z,Y)*v(Y)+R(Z,Z)*v(Z));
}
Matrix33 operator*(const Matrix33 &R, const Matrix33 &S) {
  Matrix33 T;
  forXYZ(i) forXYZ(j) {
    T(i,j)=0; forXYZ(k) T(i,j)+=R(i,k)*S(k,j);
  }
  return T;
}
float det(const Matrix33 &R) {
  return R(X,X)*(R(Z,Z)*R(Y,Y)-R(Z,Y)*R(Y,Z))
    -R(Y,X)*(R(Z,Z)*R(X,Y)-R(Z,Y)*R(X,Z))
    +R(Z,X)*(R(Y,Z)*R(X,Y)-R(Y,Y)*R(X,Z));
}
Matrix33 operator/(float c, const Matrix33 &R) {
  Matrix33 T;
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
  vector<Vector3> r; // vertices in local coordinates
  vector<vector<int> > faces;
  // properties of the body /////////////////////////
  bool locked;       // if set true, don't integrate
  float m;           // mass
  InertiaTensor I;   // moments of inertia (3x3 matrix)
  // state of the body //////////////////////////////
  Vector3 p;         // position of the center of mass
  Vector3 K;         // momentum
  Matrix33 R;        // orientation
  Vector3 L;         // angular momentum
  // auxiliary variables ///////////////////////////
  float inv_m;       // 1/m
  Matrix33 inv_I;    // 1/I
  Vector3 F;
  Vector3 tau;
  Vector3 v;         // velocity
  Vector3 omega;     // angular velocity
  vector<Vector3> vertices; // rotated and shifted r's
  // forces and constraints
  // ...
  Body() {
    m = 1.0;
    I(X,X)=I(Y,Y)=I(Z,Z)=m;
    inv_m = 1.0/m;
    inv_I = 1.0/I;
  };
  void update_vertices();
  void integrator(float dt);
  void loadObj(const string & file);
};

/**
 * rotate and shift all vertices from local to universe
 */
void Body::update_vertices() {
  vertices.resize(r.size());
  for(int i=0; i<r.size(); i++) vertices[i]=R*r[i]+p;
}

/**
 * Euler integrator
 */
void Body::integrator(float dt) {
  v     = inv_m*K;
  omega = inv_I*L;
  // p     = p + v*dt;                // shift
  K     = K + F*dt;                   // push
  R     = Rotation(omega*dt)*R;           // rotate
  // L     = L + tau*dt;                 // spin
  update_vertices();
}

/**
 * class that stores all bodies, forces and constraints
 */ 
class Universe {
public:
  list<Body*> bodies;
  // create universe by popultaing bodies
  Universe() {
    // ...
  }
  // evolve universe
  void evolve(float dt) {
    // compute force F and torque tau for each object
    for(list<Body*>::iterator body=bodies.begin();
	body!=bodies.end(); body++) {
      // ..
    }
    // integrate
    for(list<Body*>::iterator body=bodies.begin();
	body!=bodies.end(); body++) 
      if(!(*body)->locked) {
	(*body)->integrator(dt);
      }
  }
};

/**
 * Auxiliary functions translate moments of Inertia
 */
 InertiaTensor dI(float m, const Vector3 &r) {
  InertiaTensor I;
  float r2 = r*r;
  forXYZ(j) forXYZ(k) I(j,k) = m*((j==k)?r2:0-r(j)*r(k));
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
  Vector3 da = a.p-c.p;
  Vector3 db = b.p-c.p;
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
	  r.push_back(Vector3(x,y,z));
	} else if (initialVal=="f") {
	  int v1, v2, v3;
	  instream >> v1 >> v2 >> v3;
	  vector<int> triangle;
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


/*** GLUT **/
/**
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
  cout << frame << endl;
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
    Vector3 &pos = (*body)->p;
    glTranslatef(pos.v[X], pos.v[Y], pos.v[Z]);
    glBegin(GL_LINES);
    for(int i=0; i<(*body)->faces.size(); i++) {      
      for(int j=0; j<(*body)->faces[i].size(); j++) {
	int k = (*body)->faces[i][j]; 
	glVertex3fv((*body)->vertices[k].v);
      }
      int k = (*body)->faces[i][0]; 
      glVertex3fv((*body)->vertices[k].v);
    }
    glEnd();
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
  b.F = Vector3(0,-10000,0);
  b.loadObj("assets/sphere.obj");
  b.L = Vector3(0.01,0.0,0.0);
  myuniverse.bodies.push_back(&b);
  
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
