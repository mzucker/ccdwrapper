#include "CCDWrapper.h"
#include "GlCamera.h"
#include <sstream>
#include <iostream>

#ifdef __APPLE__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Glut/glut.h>
#else
#include <GL/glut.h>
#endif


using namespace ccdw;

enum { ncolors = 6 };

static const vec3 ccolors[ncolors] = {
    vec3(1.0, 0.0, 0.0),
    vec3(1.0, 1.0, 0.0),
    vec3(0.0, 1.0, 0.0),
    vec3(0.0, 1.0, 1.0),
    vec3(0.0, 0.0, 1.0),
    vec3(1.0, 0.0, 1.0),
};

//////////////////////////////////////////////////////////////////////

class CCDDemo;

CCDDemo* _instance = 0;

void demo_reshape(int w, int h);
void demo_display();
void demo_keyboard(unsigned char key, int x, int y);
void demo_mouse(int button, int state, int x, int y);
void demo_motion(int x, int y);
void demo_timer(int value);

GlCamera::MouseMode btn2cam(int button) {
  switch (button) {
  case GLUT_LEFT_BUTTON: {
    int mod = glutGetModifiers();
    switch (mod) {
    case 0:
      return GlCamera::MOUSE_ROTATE;
    case GLUT_ACTIVE_CTRL:
      return GlCamera::MOUSE_ZOOM;
    case GLUT_ACTIVE_SHIFT:
      return GlCamera::MOUSE_PAN_XY;
    default:
      break;
    }
    return GlCamera::MOUSE_NONE;
  }
  case GLUT_MIDDLE_BUTTON:
    return GlCamera::MOUSE_ZOOM;
  case GLUT_RIGHT_BUTTON:
    return GlCamera::MOUSE_PAN_XY;
  default:
    break;
  }
  return GlCamera::MOUSE_NONE;
}

class CCDDemo {
public:

    int width, height;
    
    GlCamera camera;
    GlCamera::MouseMode mmode;

    ccd_real_t arena_radius;
    std::vector<vec3> points;
    std::vector<TransformedConvex*> objects;

    std::vector<vec3> pos_rate;
    std::vector<vec3> rot_rate;

    std::vector<Report> reports;

    std::vector< std::vector<vec3> > opoints;

    std::vector< bool > colliding;

    Checker checker;
  
    QueryType qtype;
    ccd_real_t dmin;

    DrawHelper helper;

    bool animating;
    bool draw_points;
    bool draw_spheres;
    bool wireframe;

    static double randReal() {
        return rand() / double(RAND_MAX);
    }

    static vec3 randVec() {
        return vec3( randReal()*2-1,
                     randReal()*2-1,
                     randReal()*2-1 );
    }

    void setupBasicLight(const vec4& p) {

        glMatrixMode(GL_PROJECTION);
        camera.loadMatrix(GlCamera::MATRIX_PROJECTION);

        glMatrixMode(GL_MODELVIEW);
        camera.loadMatrix(GlCamera::MATRIX_MODELVIEW);

        const GLfloat lightdim[4] = { 0.3, 0.3, 0.3, 1.0 };
        const GLfloat lightbrt[4] = { 0.7, 0.7, 0.7, 1.0 };
        const GLfloat white[4] = { 1, 1, 1, 1 };
        GLfloat position[4] = { p[0], p[1], p[2], p[3] };

        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_LIGHT0);
  
        glLightfv(GL_LIGHT0, GL_POSITION, position);
        glLightfv(GL_LIGHT0, GL_AMBIENT,  lightdim);
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightbrt);
        glLightfv(GL_LIGHT0, GL_SPECULAR, white);
  
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 120.0f);
  
        glEnable(GL_LIGHTING);


    }
        

    CCDDemo(int argc, char** argv) {

        _instance = this;

        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH |
                            GLUT_RGB | GLUT_MULTISAMPLE);

        width = 640;
        height = 480;
        mmode = GlCamera::MOUSE_NONE;
        
        glutInitWindowSize(width, height);
        
        glutCreateWindow("CCD Demo");

        glutReshapeFunc(demo_reshape);
        glutDisplayFunc(demo_display);
        glutKeyboardFunc(demo_keyboard);
        glutMouseFunc(demo_mouse);
        glutMotionFunc(demo_motion);

        glClearColor(1,1,1,1);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glFrontFace(GL_CCW);
        

        camera.aim(vec3(0, 0, 6),
                   vec3(0, 0, 0),
                   vec3(0, 1, 0));

        camera.setPerspective();

        camera.setHomePosition();

        setupBasicLight(vec4(1,1,1,0));
        
        cubePoints(20, points);

        arena_radius = 1.5;

        objects.push_back(transform(dilate(new Box(vec3(0.5, 0.5, 0.5)), 0.25), Transform3(vec3(1, 0, 0))));
        objects.push_back(transform(capsule(0.5, 0.25), Transform3(vec3(0, -1, 0))));
        objects.push_back(transform(new Box(vec3(1, 1, 1)), Transform3(vec3(-1, 0, 0))));
        objects.push_back(transform(sphere(0.5), Transform3(vec3(0, 1, 0))));
        objects.push_back(transform(new Cylinder(1, 0.5), Transform3(vec3(0, 0, -1))));


        /*
          CCD_BOX(box);
          box.pos.v[0] = -1;
          box.x = box.y = box.z = 1;

          objects.push_back(new CCDConvex(box));

          box.pos.v[0] = 1;

          objects.push_back(new CCDConvex(box));
        */

        for (size_t i=0; i<objects.size(); ++i) {
            pos_rate.push_back( 0.02 * randVec() );
            rot_rate.push_back( 0.02 * randVec() );
        }

        colliding.resize( objects.size(), false );
    
        animating = false;
        draw_points = false;
        draw_spheres = false;
        wireframe = true;

        qtype = QUERY_PENETRATION;
        dmin = 3;

        checkAll();
    
        glutTimerFunc(20, demo_timer, 0);

    }

    void checkAll() {

        colliding.clear();
        colliding.resize(objects.size(), false);
        reports.clear();

        Report report;

        size_t queries = 0;
        //TimeStamp start = TimeStamp::now();

        for (size_t i=0; i<objects.size(); ++i) {
            Convex* c1 = objects[i];
            for (size_t j=0; j<i; ++j) {
                Convex* c2 = objects[j];
                std::cerr << "querying " << *c1 << " against " << *c2 << "\n";
                bool collides = checker.query(qtype, &report, c1, c2, dmin);
                ++queries;
                if (collides) { 
                    colliding[i] = colliding[j] = true;
                }
                if (report.flags & (HAVE_POSITION | HAVE_SEPARATION)) {
                    reports.push_back(report);
                }
            }
        }

        //TimeStamp end = TimeStamp::now();
        //double elapsed = (end-start).toDouble();
        double elapsed = 0;

        if (1) {
            std::cout << "did " << queries << " queries in " << elapsed << "s. "
                      << "(" << (elapsed/queries) << " per query)\n";
        }

    }

    static quat fromOmega(const vec3& v) {
        double l = v.norm();
        vec3 vn = v / l; 
        return quat(Eigen::AngleAxisd(l, vn));
    }

    virtual void timer(int value) {

        if (animating) {

            for (size_t i=0; i<objects.size(); ++i) {

                TransformedConvex* c = objects[i];

                quat q = c->xform.rotation();
                vec3 p = c->xform.translation();

                p += pos_rate[i];
                q = fromOmega(rot_rate[i]) * q;

                for (int axis=0; axis<3; ++axis) {
                    if ( (p[axis] > arena_radius && pos_rate[i][axis] > 0) ||
                         (p[axis] < -arena_radius && pos_rate[i][axis] < 0) ) {
                        pos_rate[i][axis] *= -1;
                        rot_rate[i][axis] *= -1;
                    }
                }

                c->xform = Transform3(q, p);

            }

            checkAll();

        }

        glutPostRedisplay();
        glutTimerFunc(20, demo_timer, 0);

    }

    size_t findShape(const Convex* c) { 
        for (size_t i=0; i<objects.size(); ++i) {
            if (static_cast<const Convex*>(objects[i]) == c) {
                return i;
            }
        }
        return -1;
    }

    void drawString(int rasterx, int rastery, 
                    const std::string& str,
                    void* font=GLUT_BITMAP_8_BY_13,
                    int linespacing=13) {

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(0, width, height, 0);
  
        glRasterPos2f(rasterx, rastery);
        for (std::string::const_iterator i=str.begin(); i!=str.end(); ++i) {
            if (*i == '\n') {
                rastery += linespacing;
                glRasterPos2f(rasterx, rastery);
            } else {
                glutBitmapCharacter(font, *i);
            }
        }

        glPopMatrix();

        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();

    }
    
    void viewChanged(bool interactive) {
        
        glMatrixMode(GL_MODELVIEW);
        camera.loadMatrix(GlCamera::MATRIX_MODELVIEW);
        
        glutPostRedisplay();
        
    }
        

    void reshape(int w, int h) {

        width = w;
        height = h;

        glClearColor(1,1,1,1);
  
        glViewport(0, 0, width, height);
        camera.setViewport(0, 0, width, height);

        glMatrixMode(GL_PROJECTION);
        camera.loadMatrix(GlCamera::MATRIX_PROJECTION);

        viewChanged(true);

    }

    void mouse(int button, int state, int x, int y) {
        
        if (state == GLUT_DOWN) {
            mmode = btn2cam(button);
            camera.mousePress(x,y,mmode);
        } else {
            camera.mouseRelease(x,y,mmode);
            mmode = GlCamera::MOUSE_NONE;
        }
        
        viewChanged(state == GLUT_DOWN);
    }


    void motion(int x, int y) {
        
        if (mmode != GlCamera::MOUSE_NONE) {
            camera.mouseMove(x,y,mmode);
            viewChanged(true);
        }
  
    }
    
    
    virtual void display() {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glPushAttrib(GL_POLYGON_BIT);

        if (wireframe) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        glDisable(GL_CULL_FACE);
    
        for (size_t i=0; i<objects.size(); ++i) {
            vec3 color = ccolors[i % ncolors];
            if (colliding[i]) {
                for (int i=0; i<3; ++i) {
                    if (!color[i]) { color[i] = 0.75; }
                }
                color *= 0.75;
            }
            glColor3dv(&color[0]);
            objects[i]->render(helper);
        }

        glPopAttrib();

        if (draw_spheres) {
            glColor3ub(191,191,191);
            glPushAttrib(GL_POLYGON_BIT);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            for (size_t i=0; i<objects.size(); ++i) {
                vec3 c;
                objects[i]->center(c);
                glMatrixMode(GL_MODELVIEW);
                glPushMatrix();
                glTranslated(c[0], c[1], c[2]);
                gluSphere(helper.getQuadric(), objects[i]->maxDist()+1e-3,
                          helper.slices, helper.sstacks);
                glPopMatrix();
            }
            glPopAttrib();
        }

        for (size_t i=0; i<reports.size(); ++i) {
            const Report& ri = reports[i];
            if (ri.flags & HAVE_POSITION) {
                glColor3ub(255,255,0);
            } else {
                glColor3ub(255,0,255);
            }
            glstuff::draw_cylinder(helper.getQuadric(), 
                                   ri.pos2, ri.pos1, 0.02f);
            glstuff::color(ccolors[ findShape( ri.c1 ) % ncolors ] );
            glstuff::draw_arrow(helper.getQuadric(),
                                ri.pos1, ri.pos1 - 0.3 * ri.direction, 0.04f);
            glstuff::color(ccolors[ findShape( ri.c2 ) % ncolors ] );
            glstuff::draw_arrow(helper.getQuadric(),
                                ri.pos2, ri.pos2 + 0.3 * ri.direction, 0.04f);
        }

        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);

        glColor3ub(0,0,0);
        glutWireCube(2*arena_radius);

        if (draw_points) {

            vec3 p;
            //vec3 fwd(camera.modelview().row(2).trunc());
            vec3 fwd(camera.modelview().block(2,0,1,3));

            glPointSize(2.0);
            glBegin(GL_POINTS);
            for (size_t i=0; i<objects.size(); ++i) {
                vec3 color = ccolors[ i % ncolors ] * 0.25;
                glstuff::color(color);
                const Convex* c = objects[i];
                Convex* d = NULL;
                if (dmin) {
                    c = d = dilate(c, 0.5*dmin);
                }
                for (size_t j=0; j<points.size(); ++j) {
                    if (fwd.dot(points[j]) > 0) {
                        c->support(points[j], p);
                        glstuff::vertex(p);
                    }
                }
                delete d;
            }
            glEnd();

        }
    
        glPopAttrib();

        const char* algs[] = {
            "GJK", "MPR"
        };

        const char* qtypes[] = {
            "intersect",
            "separation",
            "penetration",
        };

        std::ostringstream ostr;
        ostr << "Algorithm: " << algs[checker.algorithm] << "\n";
        ostr << "Query type: " << qtypes[qtype] << "\n";
        ostr << "Min dist: " << dmin << "\n";

        drawString(10, 20, ostr.str());

        glutSwapBuffers();

    }

    virtual void keyboard(unsigned char key, int x, int y) {
        switch (key) {
        case '\r':
        case '\n':
            animating = !animating;
            break;
        case 'p':
            draw_points = !draw_points;
            glutPostRedisplay();
            break;
        case 's':
            draw_spheres = !draw_spheres;
            glutPostRedisplay();
            break;
        case 'w':
            wireframe = !wireframe;
            glutPostRedisplay();
            break;
        case 'a':
            checker.algorithm = AlgorithmType((checker.algorithm+1)%NUM_ALGORITHMS);
            checkAll();
            glutPostRedisplay();
            break;
        case '+':
        case '=':
            dmin = std::min(dmin+ccd_real_t(0.125), 2*arena_radius);
            checkAll();
            glutPostRedisplay();
            break;
        case '-':
            dmin = std::max(dmin-ccd_real_t(0.125), ccd_real_t(0));
            checkAll();
            glutPostRedisplay();
            break;
        case 'q':
            qtype = QueryType((qtype+1)%NUM_QUERY_TYPES);
            checkAll();
            glutPostRedisplay();
            break;
        case 27: // ESC
            exit(0);
            break;
        case ' ': // SPACE
            camera.recallHomePosition(); 
            viewChanged(false);
            break;
        }
    }

    void run() {

        glutMainLoop();

    }


};

//////////////////////////////////////////////////////////////////////

void demo_reshape(int w, int h) {
    _instance->reshape(w, h);
}

void demo_display() {
    _instance->display();
}

void demo_keyboard(unsigned char key, int x, int y) {
    _instance->keyboard(key, x, y);
}

void demo_mouse(int button, int state, int x, int y) {
    _instance->mouse(button, state, x, y);
}

void demo_motion(int x, int y) {
    _instance->motion(x, y);
}

void demo_timer(int value) {
    _instance->timer(value);
}

int main(int argc, char** argv) {

    CCDDemo demo(argc, argv);
    demo.run();
  
    return 0;
}
