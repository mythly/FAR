//
//  cv_face.c
//  FAR
//
//  Created by MO TAO on 16/2/2016.
//  Copyright © 2016 MO TAO. All rights reserved.
//

#include "cv_face.h"
#include "far.h"

FART *cv_far = NULL;
int cv_count = 0;

void cv_init(const unsigned char *image_gray, int image_width, int image_height, cv_rect_t face_rect)
{
    cv_release();
    cv_far = new FART(image_gray, image_width, image_height, face_rect, &std::cout);
}

cv_rect_t cv_face_track(const unsigned char *image_gray)
{
    assert(cv_far != NULL);
    cv_rect_t rect = cv_far->track(image_gray);
    return rect;
}

bool cv_check()
{
    if (cv_far == NULL)
        return false;
    if (cv_far->error < threshold_error)
        ++cv_count;
    else
        cv_count = 0;
    if (cv_count >= 10) {
        cv_count = 0;
        return false;
    }
    return true;
}

void cv_release()
{
    if (cv_far != NULL) {
        delete cv_far;
        cv_far = NULL;
        cv_count = 0;
    }
}

