//
//  far.h
//  FAR
//
//  Created by MO TAO on 2/2/2016.
//  Copyright © 2016 MO TAO. All rights reserved.
//

#ifndef far_h
#define far_h

#include "stdlib.h"
#include "cv_common.h"

/// @brief 人脸信息结构体
typedef struct cv_face_t {
    cv_rect_t rect;			///< 代表面部的矩形区域
    float score;			///< 置信度，用于筛除负例，与人脸照片质量无关，值越高表示置信度越高。
    int yaw;			///< 水平转角，真实度量的左负右正
    int pitch;			///< 俯仰角，真实度量的上负下正
    int roll;			///< 旋转角，真实度量的左负右正
    int ID;				///< faceID，用于表示在实时人脸跟踪中的相同人脸在不同帧多次出现，在人脸检测的结果中无实际意义
} cv_face_t;

/// @brief  人脸朝向
typedef enum {
    CV_FACE_UP = 0,		///< 人脸向上，即人脸朝向正常
    CV_FACE_LEFT = 1,	///< 人脸向左，即人脸被逆时针旋转了90度
    CV_FACE_DOWN = 2,	///< 人脸向下，即人脸被逆时针旋转了180度
    CV_FACE_RIGHT = 3	///< 人脸向右，即人脸被逆时针旋转了270度
} cv_face_orientation;

cv_handle_t cv_face_create_tracker();

void cv_face_destroy_tracker(cv_handle_t tracker_handle);

cv_result_t
cv_face_track(
    cv_handle_t tracker_handle,
    const unsigned char *image,
    cv_pixel_format pixel_format,
    int image_width,
    int image_height,
    int image_stride,
    cv_face_orientation orientation,
    cv_face_t **p_faces_array,
    int *p_faces_count
);

void cv_face_reset_tracker(cv_handle_t tracker_handle);

void cv_face_release_tracker_result(cv_face_t *faces_array, int faces_count);

#endif /* far_h */
