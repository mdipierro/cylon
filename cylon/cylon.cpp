/**
 * Import the necessary libraries
 */
#include "fstream"
#include <sstream>
#include <string>
#include <iostream>
#include "math.h"
#include "vector"
#include "set"
#if defined(_MSC_VER)
#include <gl/glut.h>
#else
#include <GLUT/glut.h>
#endif
using namespace std;

/**
 * Define useful macros and constants
 */
#define array vector // to avoid name collions
#define forXYZ(i) for(int i=0; i<3; i++)
#define forEach(i,s) for(i=(s).begin();i!=(s).end();i++)
#define OBJ(iterator) (*(*iterator))
const float PRECISION = 0.00001;
const int X = 0;
const int Y = 1;
const int Z = 2;

/**
 * Define class Vector
 */
class Vector {
public:
  float v[3];
  Vector(float x=0, float y=0, float z=0) {
    v[X] = x; v[Y] = y; v[Z] = z;
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
  return Vector(v(X)-w(X),v(Y)-w(Y),v(Z)-w(Z));
}
float operator*(const Vector &v, const Vector &w) {
  return v(X)*w(X) + v(Y)*w(Y) + v(Z)*w(Z);
}
float norm(const Vector &v) {
  return sqrt(v*v);
}
Vector versor(const Vector &v) {
  float d = norm(v);
  return (d>0)?(v/d):v;
}
Vector cross(const Vector &v, const Vector &w) {
  return Vector(v(Y)*w(Z)-v(Z)*w(Y),
		v(Z)*w(X)-v(X)*w(Z),
		v(X)*w(Y)-v(Y)*w(X));
}

/**
 * Class Matrix is a base for Rotation and InertiaTensor
 */
class Matrix {
public:
  float m[3][3];
  Matrix() { forXYZ(i) forXYZ(j) m[i][j] = 0; }
  const float operator()(int i, int j) const { return m[i][j]; }
  float &operator()(int i, int j) { return m[i][j]; }
};

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
  forXYZ(i) forXYZ(j) forXYZ(k) T(i,j) += R(i,k)*S(k,j);
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
  T(X,X) = (R(Z,Z)*R(Y,Y)-R(Z,Y)*R(Y,Z))*d;
  T(X,Y) = (R(Z,Y)*R(X,Z)-R(Z,Z)*R(X,Y))*d;
  T(X,Z) = (R(Y,Z)*R(X,Y)-R(Y,Y)*R(X,Z))*d;
  T(Y,X) = (R(Z,X)*R(Y,Z)-R(Z,Z)*R(Y,X))*d;
  T(Y,Y) = (R(Z,Z)*R(X,X)-R(Z,X)*R(X,Z))*d;
  T(Y,Z) = (R(Y,X)*R(X,Z)-R(Y,Z)*R(X,X))*d;
  T(Z,X) = (R(Z,Y)*R(Y,X)-R(Z,X)*R(Y,Y))*d;
  T(Z,Y) = (R(Z,Y)*R(X,Y)-R(Z,Y)*R(X,X))*d;
  T(Z,Z) = (R(Y,Y)*R(X,X)-R(Y,X)*R(X,Y))*d;
  return T;
}

/**
 * Class Rotation is a Matrix
 */
class Rotation : public Matrix {
public:
  Rotation() { forXYZ(i) forXYZ(j) m[i][j] = (i==j)?1:0; }
  Rotation(const Vector& v) {
    float theta = norm(v);    
    if(theta<PRECISION) {
      forXYZ(i) forXYZ(j) m[i][j] = (i==j)?1:0;
    } else {
      float s = sin(theta), c = cos(theta);
      float t = 1-c;
      float x = v(X)/theta, y = v(Y)/theta, z = v(Z)/theta;
      m[X][X] = t*x*x+c;   m[X][Y] = t*x*y-s*z; m[X][Z] = t*x*z+s*y;
      m[Y][X] = t*x*y+s*z; m[Y][Y] = t*y*y+c;   m[Y][Z] = t*y*z-s*x;
      m[Z][X] = t*x*z-s*y; m[Z][Y] = t*y*z+s*x; m[Z][Z] = t*z*z+c;
    }
  }
};

/**
 * Class InertiaTensor is also a Matrix
 */
class InertiaTensor : public Matrix {};
InertiaTensor operator+(const InertiaTensor &a, const InertiaTensor &b) {
  InertiaTensor c = a;
  forXYZ(i) forXYZ(j) c(i,j) += b(i,j);
  return c;
}

/**
 * Class Body (describes a rigid body object)
 */
class Body {
public:
  // object shape
  float radius;      // all vertices inside radius;
  array<Vector> r;   // vertices in local coordinates
  array<array<int> > faces;
  Vector color;      // color of the object; 
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
    radius = 0;
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
  void loadObj(const string & file, float scale);
  void draw();
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
 * The constructor can be specific of the force.
 * The apply methods adds the contribution to F and tau
 */
class Force {
public:
  virtual void apply(float dt)=0;
  virtual void draw() {};
};

/**
 * Gravity is a Force
 */
class GravityForce : public Force {
public:  
  Body *body;
  float g;
  GravityForce(Body *body, float g = 0.01) {
    this->body = body; this->g = g;    
  }
  void apply(float dt) {
    body->F(Y) -= (body->m)*g;
  }
};

/**
 * Spring is a Force
 */
class SpringForce : public Force {
public:
  Body *bodyA;
  Body *bodyB;
  int iA, iB;
  float kappa, L;
  SpringForce(Body *bodyA, int iA, Body *bodyB, int iB,
	      float kappa, float L) {
    this->bodyA = bodyA; this->iA = iA;
    this->bodyB = bodyB; this->iB = iB;
    this->kappa = kappa; this->L = L;
  }
  void apply(float dt) {
    Vector d = bodyB->vertices[iB]-bodyA->vertices[iA];
    float n = norm(d);
    if(n>PRECISION) {
      Vector F = kappa*(n-L)*(d/n);
      bodyA->F = bodyA->F+F;
      bodyB->F = bodyB->F-F;
      bodyA->tau = bodyA->tau + cross(bodyA->Rr[iA],F);
      bodyB->tau = bodyB->tau - cross(bodyB->Rr[iB],F);
    }
  }
  void draw();
};

/**
 * A Spring can be anchored to a pin.
 */
class AnchoredSpringForce : public Force {
public:
  Body *body;
  int i;
  Vector pin;
  float kappa, L;
  AnchoredSpringForce(Body *body, int i, Vector pin,
		      float kappa, float L) {
    this->body = body; this->i = i;
    this->pin = pin;
    this->kappa = kappa; this->L = L;
  }
  void apply(float dt) {
    Vector d = body->vertices[i]-pin;
    float n = norm(d);
    Vector F = kappa*(n-L)*d/n;
    body->F = body->F+F;
    body->tau = body->tau + cross(body->Rr[i],F);
  }
};


/**
 * Friction is alos a Force
 * (this ignores the shape of the body, assumes a sphere)
 */
class FrictionForce: public Force {
public:
  Body *body;
  float gamma;
    FrictionForce(Body *body, float gamma) {
      this->body = body; this->gamma = gamma;
    }
    void apply(float dt) {
      body->F = body->F-gamma*(body->v)*dt; // ignores shape
    }
};


/**
 * class Water, one instance floods the universe
 */
class Water: public Force {
public:
  float level, wave, speed;
  float m[41][41];  
  float t, x0,dx;
  set<Body*> *bodies;
  Water(set<Body*> *bodies, float level, 
	float wave=0.2, float speed=0.2) {
    this->bodies = bodies;
    this->level = level; this->wave = wave, this->speed = speed;
    t = 0; x0 = 5.0; dx = 0.25;
  }
  void apply(float dt) {
    set<Body*>::iterator ibody;
    t = t+dt;
    for(int i=0; i<41; i++)
      for(int j=0; j<41; j++)
	this->m[i][j] = level+wave/2*sin(speed*t+i)+wave/2*sin(0.5*i*j);
    forEach(ibody,*bodies) {
      Body &body = OBJ(ibody); // dereference
      int i = (body.p(X)+x0)/dx;
      int j = (body.p(Z)+x0)/dx;
      if(body.p(Y)<m[i][j]) {
	body.F = body.F + (Vector(0,1,0)-2.0*(body.v))*dt;    
	body.L = (1.0-dt)*body.L;
      }
    }
  }
  void draw();
};

/**
 * A Constraint has detect and resolve methods.
 */
class Constraint {
public:  
  virtual bool detect()=0;
  virtual void resolve(float dt)=0;
  virtual void draw() {}
  Vector impluse(const Body &A, const Body &B, 
		 Vector &r_A, Vector &r_B, Vector n, float c) {
    float IA; // FIX
    float IB; // FIX
    Vector r_cA = A.R*r_A;
    Vector r_cB = B.R*r_B;
    Vector v_cA = cross(A.omega,r_A)+A.v;
    Vector v_cB = cross(B.omega,r_B)+B.v;
    Vector crossA = cross(r_cA,n);
    Vector crossB = cross(r_cB,n);
    Vector dF = (-(c-1)/(A.inv_m+B.inv_m+
			 crossA*crossA/IA+
			 crossB*crossB/IB)*(v_cB-v_cB)*n)*n;
  }
};

/**
 * Class to deal with collision with static plane
 * (ignore rotation)
 */
class PlaneConstraint: public Constraint {
public:
  Body *body;
  Vector n,d;
  float penetration, restitution;
  // d is the discance of the plane from origin
  // n is a versor orthogonal to plane 
  // (in opposite direction from collision)
  PlaneConstraint(Body *body, float restitution,
		  const Vector &d, const Vector &n) {
    this->body = body;
    this->n = n; this->d = d;
    this->restitution = restitution;
  }
  bool detect() {
    penetration = body->p*n-d*n + body->radius;
    return penetration>=0;
  }
  void resolve(float dt) {
    // move the object back is stuck on plane 
    body->p = body->p - penetration*n;
    float K_ortho = n*body->K;
    // optional, deal with friction 
    Vector L_ortho = -(body->radius)*cross(n,body->K-K_ortho);
    body->L = (n*body->L)*n + L_ortho;
    // reverse momentum
    if(K_ortho>0)
      body->K = body->K - (restitution+1)*(K_ortho)*n;
  }
};

/**
 * Default all-2-all collision handler
 */
class All2AllCollisions : public Constraint {
public:  
  float restitution;
  set<Body*> *bodies;
  All2AllCollisions(set<Body*> *bodies, float c=0.5) {
    this->bodies=bodies;
    restitution = c;
  }
  bool detect() { return true; }
  void All2AllCollisions::resolve(float dt) {
    set<Body*>::iterator ibodyA, ibodyB;
    forEach(ibodyA,*bodies) {
      forEach(ibodyB,*bodies) {
	Body &A = OBJ(ibodyA); // dereference
	Body &B = OBJ(ibodyB); // dereference
	Vector d = B.p-A.p;
	Vector v_closing = B.v-A.v;
	float penetration = A.radius+B.radius-norm(d);	  
	if(penetration>0 && v_closing*d<0) {
	  Vector q = (penetration/(A.m+B.m))*versor(d);
	  A.p = A.p-B.m*q;
	  B.p = B.p+A.m*q;
	  Vector impulse = (-restitution*A.m*B.m/(A.m+B.m))*v_closing;
	  A.K = A.K-impulse;
	  B.K = B.K+impulse;
	}
      } 
    }
  }
};


/**
 * A Universe stores bodies, forces, constraints
 * and evolves in time.
 */ 
class Universe {
public:
  float dt;
  // universe state
  set<Body*> bodies;
  set<Force*> forces;
  set<Constraint*> constraints;
  // useful iterators
  set<Body*>::iterator ibody;
  set<Force*>::iterator iforce;
  set<Constraint*>::iterator iconstraint;
  int frame;
  Universe() {
    frame = 0;
  }
  ~Universe() { 
    forEach(ibody,bodies) delete (*ibody);
    forEach(iforce,forces) delete (*iforce);
    forEach(iconstraint,constraints) delete (*iconstraint);
  }
  // evolve universe
  void evolve() {    
    // clear forces and troques
    forEach(ibody,bodies)
      OBJ(ibody).clear();
    // compute forces and torques
    forEach(iforce,forces)
      OBJ(iforce).apply(dt);
    callback();
    // integrate
    forEach(ibody,bodies)
      if(!OBJ(ibody).locked)
	OBJ(ibody).integrator(dt);
    // handle collisions (not quite right yet)
    forEach(iconstraint,constraints)
      if(OBJ(iconstraint).detect())
	OBJ(iconstraint).resolve(dt);
    frame++;
  }
public:
  virtual void build_universe()=0;
  virtual void callback()=0;
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
  c.L = (a.p-c.p)*a.K+(b.p-c.p)*b.K;
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

/**
 * Function that loads an wavefront obj file into a Body.
 */
void Body::loadObj(const string & file,float scale=0.5) {
  ifstream input;
  string line;
  float x,y,z;
  input.open(file.c_str());
  if(input.is_open()) {
    while(input.good()) {
      getline(input, line);
      if(line.length()>0) {
	string initialVal;
	istringstream instream;
	instream.str(line);
	instream >> initialVal;
	if(initialVal=="v") {
	  instream >> x >> y >> z;
	  Vector p = scale*Vector(x,y,z);
	  r.push_back(p);
	  radius = max(radius,norm(p));
	} else if (initialVal=="f") {
	  array<int> path;
	  while(instream >> x) path.push_back(x-1);
	  faces.push_back(path);
	}
      }
    }
    update_vertices();
  }
}

/**
 * Make My Universe!
 */
class MyUniverse : public Universe {
public:
  void build_universe() {
    Body *b_old = 0;
    for(int i=0; i<4; i++) {
      Body &b = *new Body();
      b.color=Vector(((i+1)%4)?1:0,i%2,(i%3)?1:0);
      b.loadObj("assets/sphere.obj");
      b.p = Vector(i,i+2,-i);
      b.K = Vector(0.1*i,0.01*i,0);
      b.L = Vector(0.5,0.5*i,0.1*i);
      bodies.insert(&b);
      forces.insert(new GravityForce(&b,0.01));
      constraints.insert(new PlaneConstraint(&b,0.9,Vector(0,0,0),
					     Vector(0,-1,0)));
      constraints.insert(new PlaneConstraint(&b,0.9,Vector(5,0,0),
					     Vector(1,0,0)));
      constraints.insert(new PlaneConstraint(&b,0.9,Vector(-5,0,0),
					     Vector(-1,0,0)));
      constraints.insert(new PlaneConstraint(&b,0.9,Vector(0,0,5),
					     Vector(0,0,1)));
      constraints.insert(new PlaneConstraint(&b,0.9,Vector(0,0,-5),
					     Vector(0,0,-1)));
      // forces.insert(new FrictionForce(&b,0.5));
      if(i==2)
	forces.insert(new SpringForce(b_old,0,&b,0,0.01,0));
      b_old = &b;
    }
    constraints.insert(new All2AllCollisions(&bodies,0.8));
    // forces.insert(new Water(&bodies,3.0));
  }
  void callback() {}
};

MyUniverse universe;

/**
 * GLUT code below
 * Creates a window in which to display the scene.
 */
void createWindow(const char* title) {
  int width = 640, height = 480;
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

/**
 * Called each frame to update the 3D scene. Delegates to
 * the application.
 */
void update() {
  // evolve world
  universe.dt = 0.016f;	// 60fps fixed rate.
  universe.evolve();
  glutPostRedisplay();
}

/**
 * Function called each frame to display the 3D scene. 
 * It draws all bodies, forces and constraints.
 */
void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(0.0,3.5,10.0, 0.0,3.5,0.0, 0.0,1.0,0.0);
  forEach(universe.ibody,universe.bodies)
    OBJ(universe.ibody).draw();
  forEach(universe.iforce,universe.forces)
    OBJ(universe.iforce).draw();
  forEach(universe.iconstraint,universe.constraints)
    OBJ(universe.iconstraint).draw();
  // update the displayed content
  glFlush();
  glutSwapBuffers();
}

/**
 * Code that draws an Body
 */
void Body::draw() {
  glPolygonMode(GL_FRONT,GL_FILL);
  for(int i=0; i<faces.size(); i++) {      
    float k = 0.5*(1.0+(float)(i+1)/faces.size());
    glColor3f(color(0)*k,color(1)*k,color(2)*k);    
    glBegin(GL_POLYGON);        
    for(int j=0; j<faces[i].size(); j++)
      glVertex3fv(vertices[faces[i][j]].v);
    glVertex3fv(vertices[faces[i][0]].v);
    glEnd();
  }
}

/**
 * Code that draws a Spring
 */
void SpringForce::draw() {
  glColor3f(0,0,0);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  glBegin(GL_LINES);
  glVertex3fv(bodyA->vertices[iA].v);
  glVertex3fv(bodyB->vertices[iB].v);
  glEnd();
}

/**
 * Draw the the Water
 */
void Water::draw() {
  glPolygonMode(GL_FRONT,GL_FILL);
  for(int i=0; i<40; i++) {
    glBegin(GL_POLYGON);
    for(int j=0; j<40; j++) {
      glColor3f(0,0,0.5+0.25*(m[i][j]-level+wave)/wave);
      glVertex3fv(Vector(-x0+dx*i,m[i][j],-x0+dx*j).v);
      glVertex3fv(Vector(-x0+dx*i,m[i][j+1],-x0+dx*j+dx).v);
      glVertex3fv(Vector(-x0+dx*i+dx,m[i+1][j+1],-x0+dx*j+dx).v);
      glVertex3fv(Vector(-x0+dx*i+dx,m[i+1][j],-x0+dx*j).v);
    }
    glEnd();
  }      
}

/**
 * Function called when the display window changes size.
 */
void reshape(int width, int height) {
  glViewport(0, 0, width, height);
}

/**
 * Function called when a mouse button is pressed.
 */
void mouse(int button, int state, int x, int y) { }

/**
 * Function called when a key is pressed.
 */
void keyboard(unsigned char key, int x, int y) {
  // '1' kick ball 1, '2' kicks ball 2, etc.
  int i = 0;
  set<Body*>::iterator ibody;
  forEach(ibody,universe.bodies) {
    if(key-49==i) {
      Body &body = OBJ(ibody); // dereference
      body.K = Vector(0.2,0.2,0);
      body.L = Vector(0,0,0.04);
    }
    i++;
  }
}

/**
 * Called when the mouse is dragged.
 */
void motion(int x, int y) { }

/**
 * The main function. Everythign starts here.
 */
int main(int argc, char** argv) {
  // Create the application and its window
  glutInit(&argc, argv);
  createWindow("Cylon");
  // fill universe with stuff
  universe.build_universe();
  // Set up the appropriate handler functions
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutDisplayFunc(display);
  glutIdleFunc(update);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);  
  // Enter into a loop
  glutMainLoop();  
  return 0;
}
