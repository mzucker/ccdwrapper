#ifndef _GLSTUFF_H_
#define _GLSTUFF_H_

#ifdef __APPLE__
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef HAVE_GLEW
#include <GL/glew.h>
#else
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#endif

#include "geom.h"

namespace glstuff {

    void check_opengl_errors(const char* context);

    size_t pwr2(size_t x);

    inline void get(GLenum which, GLdouble* data) {
        glGetDoublev(which, data);
    }

    inline void get(GLenum which, GLfloat* data) {
        glGetFloatv(which, data);
    }

    inline void translate(GLfloat x, GLfloat y, GLfloat z) {
        glTranslatef(x,y,z);
    }

    inline void translate(GLdouble x, GLdouble y, GLdouble z) {
        glTranslated(x,y,z);
    }

    inline void translate(const vec3& v) { glTranslated( v[0], v[1], v[2] ); }

    inline void color(const vec3& v) { glColor3dv( v.data() ); }
    inline void color(const vec4& v) { glColor4dv( v.data() ); }

    inline void normal(const vec3& v) { glNormal3dv( v.data() ); }

    inline void vertex(const vec3& v) { glVertex3dv( v.data() ); }

    inline void mat4_to_gl(const mat4& mat, GLdouble data[16]) {

        data[0]  = mat(0,0);
        data[1]  = mat(1,0);
        data[2]  = mat(2,0);
        data[3]  = mat(3,0);

        data[4]  = mat(0,1);
        data[5]  = mat(1,1);
        data[6]  = mat(2,1);
        data[7]  = mat(3,1);

        data[8]  = mat(0,2);
        data[9]  = mat(1,2);
        data[10] = mat(2,2);
        data[11] = mat(3,2);

        data[12] = mat(0,3);
        data[13] = mat(1,3);
        data[14] = mat(2,3);
        data[15] = mat(3,3);

    }
  

    inline void gl_to_mat4(const GLdouble data[16], mat4& mat) {

        mat(0,0) = data[0];
        mat(1,0) = data[1];
        mat(2,0) = data[2];
        mat(3,0) = data[3];

        mat(0,1) = data[4];
        mat(1,1) = data[5];
        mat(2,1) = data[6];
        mat(3,1) = data[7];

        mat(0,2) = data[8];
        mat(1,2) = data[9];
        mat(2,2) = data[10];
        mat(3,2) = data[11];

        mat(0,3) = data[12];
        mat(1,3) = data[13];
        mat(2,3) = data[14];
        mat(3,3) = data[15];

    }

    inline void load_matrix(GLfloat* data)  { glLoadMatrixf(data); }
    inline void load_matrix(GLdouble* data) { glLoadMatrixd(data); }

    inline void mult_matrix(GLfloat* data)  { glMultMatrixf(data); }
    inline void mult_matrix(GLdouble* data) { glMultMatrixd(data); }
  
    inline void get_mat4(GLenum which, mat4& mat) {
        GLdouble data[16];
        get(which, data);
        gl_to_mat4(data, mat);
    }

    inline void load_mat4(const mat4& mat) {
        GLdouble data[16];
        mat4_to_gl(mat, data);
        load_matrix(data);
    }

    inline void mult_mat4(const mat4& mat) {
        GLdouble data[16];
        mat4_to_gl(mat, data);
        mult_matrix(data);
    }

    inline void load_transform(const Transform3& tx) {
        load_mat4(tx.matrix());
    }

    inline void mult_transform(const Transform3& tx) {
        mult_mat4(tx.matrix());
    }

    enum { 
        DEFAULT_SLICES = 32,
        DEFAULT_CSTACKS = 4,
        DEFAULT_SSTACKS = 24
    };

    inline bool setup_cylinder(const vec3& p0,
                               const vec3& p1,
                               double& len) {
    
        vec3 z = p1-p0;
        len = z.norm();

        if (len <= 1e-9) {
            translate(p0);
            return false;
        }

        quat q = quatFromOneVector(z);

        mult_transform(Transform3(q, p0));

        return true;

    }
    
    inline vec3 sph2crt(double lat, double lon) {
        return vec3( sin(lat)*cos(lon),
                     sin(lat)*sin(lon),
                     cos(lat) );
    }


    inline void draw_round_box(const Box3& box,
                               double r,
                               int slices, int stacks) {

        double x0 = box.p0.x(), y0 = box.p0.y(), z0 = box.p0.z();
        double x1 = box.p1.x(), y1 = box.p1.y(), z1 = box.p1.z();

        typedef vec3 vec3;

        vec3 normals[6] = {
            vec3(-1,  0,  0),
            vec3( 1,  0,  0),
            vec3( 0, -1,  0),
            vec3( 0,  1,  0),
            vec3( 0,  0, -1),
            vec3( 0,  0,  1)
        };

        vec3 quads[6][4] = {
            { vec3(x0, y0, z0), vec3(x0, y0, z1), vec3(x0, y1, z1), vec3(x0, y1, z0) },
            { vec3(x1, y0, z0), vec3(x1, y1, z0), vec3(x1, y1, z1), vec3(x1, y0, z1) },
            { vec3(x0, y0, z0), vec3(x1, y0, z0), vec3(x1, y0, z1), vec3(x0, y0, z1) },
            { vec3(x0, y1, z0), vec3(x0, y1, z1), vec3(x1, y1, z1), vec3(x1, y1, z0) },
            { vec3(x0, y0, z0), vec3(x0, y1, z0), vec3(x1, y1, z0), vec3(x1, y0, z0) },
            { vec3(x0, y0, z1), vec3(x1, y0, z1), vec3(x1, y1, z1), vec3(x0, y1, z1) },
        };


        glBegin(GL_TRIANGLES);
        glEdgeFlag(GL_TRUE);

        for (int i=0; i<6; ++i) {

            normal(normals[i]);

            vertex(quads[i][0] + normals[i]*r);
            glEdgeFlag(GL_FALSE);
            vertex(quads[i][1] + normals[i]*r);
            glEdgeFlag(GL_TRUE);
            vertex(quads[i][3] + normals[i]*r);

            vertex(quads[i][1] + normals[i]*r);
            vertex(quads[i][2] + normals[i]*r);
            glEdgeFlag(GL_FALSE);
            vertex(quads[i][3] + normals[i]*r);
            glEdgeFlag(GL_TRUE);
            
        }


        int nst = stacks/2;
        int nsl = slices/4;

        vec3 corners[4] = {
            vec3(x1, y1, z1),
            vec3(x0, y1, z1),
            vec3(x0, y0, z1),
            vec3(x1, y0, z1),
        };

        for (int rot=0; rot<4; ++rot) {
      
            double lr = M_PI/2 * rot;

            vec3 c1 = corners[rot];
            vec3 c0 = corners[(rot+3)%4];
            vec3 cz(c1.x(), c1.y(), z0);

            for (int sl=0; sl<nsl; ++sl) {
        
                double lo0 = lr + (sl*2*M_PI)/slices;
                double lo1 = lr + ((sl+1 == nsl) ? M_PI/2 : ((sl+1)*2*M_PI)/slices);
        
                vec3 s0 = sph2crt((M_PI/2), lo0);
                vec3 s1 = sph2crt((M_PI/2), lo1);
        
                normal(s1);
                vertex(c1 + r*s1);
                normal(s0);
                glEdgeFlag(GL_FALSE);
                vertex(c1 + r*s0);
                glEdgeFlag(GL_TRUE);
                normal(s1);
                vertex(cz + r*s1);
        
                normal(s0);
                vertex(c1 + r*s0);
                normal(s0);
                vertex(cz + r*s0);
                normal(s1);
                glEdgeFlag(GL_FALSE);
                vertex(cz + r*s1);
                glEdgeFlag(GL_TRUE);

            }

            for (int bot=0; bot<2; ++bot) {

                vec3 cb1 = c1;
                vec3 cb0 = c0;

                if (bot) { 
                    cb1.z() = cb0.z() = z0;
                }

                double lb = bot*M_PI/2;
        
                for (int st=0; st<nst; ++st) {

                    double la0 = lb + (st*M_PI)/stacks;
                    double la1 = lb + ((st+1 == nst) ? M_PI/2 : ((st+1)*M_PI)/stacks);

                    for (int sl=0; sl<nsl; ++sl) {

                        int sll = sl;

                        double lo0 = lr + (sll*2*M_PI)/slices;
                        double lo1 = lr + ((sll+1 == nsl) ? M_PI/2 : ((sll+1)*2*M_PI)/slices);

                        vec3 s00 = sph2crt(la0, lo0);
                        vec3 s10 = sph2crt(la1, lo0);
                        vec3 s01 = sph2crt(la0, lo1);
                        vec3 s11 = sph2crt(la1, lo1);

                        normal(s01);
                        vertex(cb1 + r*s01);
                        normal(s00);
                        glEdgeFlag(GL_FALSE);
                        vertex(cb1 + r*s00);
                        glEdgeFlag(GL_TRUE);
                        normal(s11);
                        vertex(cb1 + r*s11);

                        normal(s00);
                        vertex(cb1 + r*s00);
                        normal(s10);
                        vertex(cb1 + r*s10);
                        normal(s11);
                        glEdgeFlag(GL_FALSE);
                        vertex(cb1 + r*s11);
                        glEdgeFlag(GL_TRUE);

                    }
      
                    vec3 s0 = sph2crt(la0, lr);
                    vec3 s1 = sph2crt(la1, lr);
        
                    normal(s1);
                    vertex(cb1 + r*s1);
                    normal(s0);
                    glEdgeFlag(GL_FALSE);
                    vertex(cb1 + r*s0);
                    glEdgeFlag(GL_TRUE);
                    normal(s1);
                    vertex(cb0 + r*s1);

                    normal(s0);
                    vertex(cb1 + r*s0);
                    normal(s0);
                    vertex(cb0 + r*s0);
                    normal(s1);
                    glEdgeFlag(GL_FALSE);
                    vertex(cb0 + r*s1);
                    glEdgeFlag(GL_TRUE);

      
                }

            }

        }

        glEnd();

    }


    inline void draw_cylinder(GLUquadric* q,
                              const vec3& p0,
                              const vec3& p1, 
                              double rad,
                              int slices=DEFAULT_SLICES,
                              int cstacks=DEFAULT_CSTACKS,
                              int cloops=1) {

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        double len;

        if (!setup_cylinder(p0, p1, len)) {
    
            return;
    
        } else {
    
            gluCylinder(q, rad, rad, len, slices, cstacks);
            gluQuadricOrientation(q, GLU_INSIDE);
            gluDisk(q, 0, rad, slices, cloops);

            translate(0, 0, len);
            gluQuadricOrientation(q, GLU_OUTSIDE);
            gluDisk(q, 0, rad, slices, cloops);


        }

        glPopMatrix();


    }

    inline void draw_box(const Box3& box) {

        double x0 = box.p0.x(), y0 = box.p0.y(), z0 = box.p0.z();
        double x1 = box.p1.x(), y1 = box.p1.y(), z1 = box.p1.z();

        typedef vec3 vec3;

        vec3 normals[6] = {
            vec3(-1,  0,  0),
            vec3( 1,  0,  0),
            vec3( 0, -1,  0),
            vec3( 0,  1,  0),
            vec3( 0,  0, -1),
            vec3( 0,  0,  1)
        };

        vec3 quads[6][4] = {
            { vec3(x0, y0, z0), vec3(x0, y0, z1), vec3(x0, y1, z1), vec3(x0, y1, z0) },
            { vec3(x1, y0, z0), vec3(x1, y1, z0), vec3(x1, y1, z1), vec3(x1, y0, z1) },
            { vec3(x0, y0, z0), vec3(x1, y0, z0), vec3(x1, y0, z1), vec3(x0, y0, z1) },
            { vec3(x0, y1, z0), vec3(x0, y1, z1), vec3(x1, y1, z1), vec3(x1, y1, z0) },
            { vec3(x0, y0, z0), vec3(x0, y1, z0), vec3(x1, y1, z0), vec3(x1, y0, z0) },
            { vec3(x0, y0, z1), vec3(x1, y0, z1), vec3(x1, y1, z1), vec3(x0, y1, z1) },
        };

        glBegin(GL_QUADS);
        for (int i=0; i<6; ++i) {
            normal(normals[i]);
            for (int j=0; j<4; ++j) {
                vertex(quads[i][j]);
            }
        }
        glEnd();

    }

    inline void draw_capsule(GLUquadric* q,
                             const vec3& p0,
                             const vec3& p1,
                             double rad,
                             int slices=DEFAULT_SLICES,
                             int sstacks=DEFAULT_SSTACKS,
                             int cstacks=DEFAULT_CSTACKS,
                             int clipplane=GL_CLIP_PLANE0) {

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        double len;

        bool doclip = (clipplane >= GL_CLIP_PLANE0);

        if (!setup_cylinder(p0, p1, len)) {
    

            gluSphere(q, rad, slices, sstacks);

    
        } else {

            if (doclip) {
                vec4 cp(0, 0, -1, 0);
                glEnable(clipplane);
                glClipPlane(clipplane, &cp[0]);
            }
    
            gluSphere(q, rad, slices, sstacks);

            if (doclip) {
                glDisable(clipplane);
            }
          

            gluCylinder(q, rad, rad, len, slices, cstacks);
            translate(0, 0, len);

            if (doclip) {
                vec4 cp(0, 0, 1, 0);
                glEnable(clipplane);
                glClipPlane(clipplane, &cp[0]);
            }

            gluSphere(q, rad, slices, sstacks);

            if (doclip) {
                glDisable(clipplane);
            }

        }

        glPopMatrix();
  
    }

    inline void draw_arrow(GLUquadric* q,
                           const vec3& p0,
                           const vec3& p1,
                           double cyl_rad,
                           double head_rad=0, // defaults to 2*cyl_rad
                           double head_len=0, // defaults to 2*cyl_rad
                           int slices=DEFAULT_SLICES,
                           int cstacks=DEFAULT_CSTACKS,
                           int hstacks=1,
                           int hloops=1,
                           int cloops=1) {

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        double len;
        if (setup_cylinder(p0, p1, len)) { 

            if (!head_rad) { head_rad = 2*cyl_rad; }
            if (!head_len) { head_len = 4*cyl_rad; }

            double clen = len - 0.5 * head_len;
            double flange_inner = 0;

            if (clen > 0) {

                // cylinder
                gluCylinder(q, cyl_rad, cyl_rad, clen, slices, cstacks);

                // base disk
                gluQuadricOrientation(q, GLU_INSIDE);
                gluDisk(q, 0, cyl_rad, slices, cloops);
      
                flange_inner = cyl_rad;

            }

            translate(0, 0, clen);

            // flange
            gluDisk(q, flange_inner, head_rad, slices, hloops);
            gluQuadricOrientation(q, GLU_OUTSIDE);

            // cone
            gluCylinder(q, head_rad, 0, head_len, slices, hstacks);


        }
  
        glPopMatrix();
  
    }

}

#endif
