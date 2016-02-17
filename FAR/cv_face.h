//
//  cv_face.h
//  FAR
//
//  Created by MO TAO on 16/2/2016.
//  Copyright © 2016 MO TAO. All rights reserved.
//

#ifndef cv_face_h
#define cv_face_h

typedef struct cv_rect_t {
    float x;
    float y;
    float width;
    float height;
} cv_rect_t;

#ifdef __cplusplus
extern "C" {
#endif
    
    void cv_init(const unsigned char *image_gray, int image_width, int image_height, cv_rect_t face_rect);
    cv_rect_t cv_face_track(const unsigned char *image_gray);
    bool cv_check();
    void cv_release();
    
#ifdef __cplusplus
}
#endif

#endif /* cv_face_h */
