#ifndef _GLCAMERA_H_
#define _GLCAMERA_H_

#include "glstuff.h"

class GlCamera {
public:

    enum CameraType {
        CAMERA_ORTHO,
        CAMERA_PERSPECTIVE
    };

    enum RotateType {
        ROTATE_2_AXIS,
        ROTATE_TRACKBALL
    };

    enum MatrixType {
        MATRIX_PROJECTION,
        MATRIX_MODELVIEW,
        MATRIX_PROJECTION_INVERSE,
        MATRIX_MODELVIEW_INVERSE,
    };

    enum MouseMode {
        MOUSE_NONE,
        MOUSE_PAN_XY,
        MOUSE_PAN_Z,
        MOUSE_ROTATE,
        MOUSE_ZOOM
    };

    enum MouseFlags {
        MOUSE_ALLOW_PAN    = 0x01,
        MOUSE_ALLOW_ROTATE = 0x02,
        MOUSE_ALLOW_ZOOM   = 0x04,
        MOUSE_ALLOW_ALL    = 0x07,
        MOUSE_ALLOW_NONE   = 0x00
    };

    //////////////////////////////////////////////////

    class PerspectiveInfo {
    public:
        float fovy;
        float near;
        float far;

        PerspectiveInfo();
        PerspectiveInfo(float fovy, float near, float far);

    };

    //////////////////////////////////////////////////

    class AimInfo {
    public:

        vec3 position;
        vec3 target;
        vec3 up;

        AimInfo();
        AimInfo(const vec3& p, const vec3& t, const vec3& u);

        // translates in the current frame
        void translate(const vec3& t);

        // rotates around the target position
        void rotate(const quat& q);

        Transform3 xform() const;

    };

    //////////////////////////////////////////////////

    GlCamera();

    //////////////////////////////////////////////////

    void setViewport(int x, int y, int w, int h);
  
    void setOrtho(const vec3& viewDims);

    void setPerspective(const PerspectiveInfo& pi);

    void setPerspective(float fovy=45.0f, 
                        float near=0.1,
                        float far=100.0);


    void setCameraType(CameraType t);
    void setRotateType(RotateType t);
  
    void aim(const AimInfo& aim);

    void aim(const vec3& position,
             const vec3& target,
             const vec3& up);

    void setZoomRange(float minSize, float maxSize);
    void setXRotRange(float xmin, float xmax);
    void setYRotRange(float ymin, float ymax);
    void set2AxisAngles(float xr, float yr);
    void setMouseFlags(int flags);

    //////////////////////////////////////////////////

    CameraType             cameraType() const;
    RotateType             rotateType() const;
    const vec3&           orthoViewDims() const;
    const PerspectiveInfo& perspectiveInfo() const;
    const AimInfo&         aimInfo() const;
    void                   get2AxisAngles(float& xr, float& yr) const;
    int                    mouseFlags() const;

    //////////////////////////////////////////////////

    void   loadMatrix(MatrixType m) const;
    void   multMatrix(MatrixType m) const;
    const mat4& getMatrix(MatrixType m) const;
  
    const Transform3& xform() const;

    const mat4& projection() const;
    const mat4& modelview() const;

    const mat4& projectionInverse() const;
    const mat4& modelviewInverse() const;

    const int* viewport() const;

    //////////////////////////////////////////////////

    void setHomePosition();
    void recallHomePosition();

    //////////////////////////////////////////////////

    void pan(float px, float py, float pz=0, bool autoscale=false);
    void pan(const vec3& p, bool autoscale=false);
    void zoom(float factor);

    //////////////////////////////////////////////////

    void mousePress(int x, int y, MouseMode m);
    void mouseMove(int x, int y, MouseMode m);
    void mouseRelease(int x, int y, MouseMode m);
    void mouseReset();

    //////////////////////////////////////////////////

    vec3 unproject(int x, int y, float z=0.5) const;
    vec3 unproject(const vec3& wincoords) const;

private:

    vec3 _hemiCoords(int x, int y) const;
    void _updateProjection();
    void _updateModelview();

    CameraType _cameraType;
    RotateType  _rotateType;

    int   _viewport[4];
    float _aspect; // derived from _viewport

    // rotating moves position but not target for trackball, just changes xrot/yrot for tilt
    // panning moves position AND target
    // zooming moves position only AND changes view dims

    class PoseInfo {
    public:
        AimInfo aim;
        float   xrot, yrot;
        vec3   orthoViewDims;
        PoseInfo();
    };

    PerspectiveInfo _pinfo;

    
    PoseInfo _home;
    PoseInfo _current;
  
    mat4       _projection; // derived from _current, _pinfo, _viewport
    float       _winsz; // derived from orthoViewDims or pinfo + aim
    float       _minsz;
    float       _maxsz;
    float       _minxrot;
    float       _maxxrot;
    float       _minyrot;
    float       _maxyrot;

    Transform3 _cxform; // derived from _current
    mat4       _modelview; // derived from _current

    mat4       _projection_inv;
    mat4       _modelview_inv;
  

    int _mouseFlags;
    bool _mouseDown;
    int  _mouseMode;
    PoseInfo _mouseOrigPose;
    int _mouseX, _mouseY;
    float _mouseWinSz;
    vec3 _mouseHemi;
  
    
};

#endif
