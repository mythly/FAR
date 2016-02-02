//
//  far.c
//  FAR
//
//  Created by MO TAO on 2/2/2016.
//  Copyright © 2016 MO TAO. All rights reserved.
//

#include "far.h"

cv_handle_t cv_face_create_tracker()
{
    return 0;
}

void cv_face_destroy_tracker(cv_handle_t tracker_handle)
{
    
}

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
              )
{
    *p_faces_array = malloc(sizeof(cv_face_t));
    cv_rect_t rand_rect;
    rand_rect.top = rand() % (image_height - image_height / 4);
    rand_rect.bottom = rand_rect.top + image_height / 4;
    rand_rect.left = rand() % (image_width - image_width / 4);
    rand_rect.right = rand_rect.left + image_width / 4;
    (*p_faces_array)[0].rect = rand_rect;
    *p_faces_count = 1;
    
    return CV_OK;
}

void cv_face_reset_tracker(cv_handle_t tracker_handle)
{
    
}

void cv_face_release_tracker_result(cv_face_t *faces_array, int faces_count)
{
    free(faces_array);
}
