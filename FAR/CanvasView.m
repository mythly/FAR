//
//  CanvasView.m
//  Created by sluin on 15/7/1.
//  Copyright (c) 2015年 SunLin. All rights reserved.
//

#import "CanvasView.h"

@implementation CanvasView
{
    CGContextRef context ;
}

- (void)drawRect:(CGRect)rect {
    [self drawFace:CGRectFromString(self.strFace) leftEye:CGPointFromString(self.strLeftEye) rightEye:CGPointFromString(self.strRightEye) mouse:CGPointFromString(self.strMouse)] ;
}

-(void)drawFace:(CGRect)face leftEye:(CGPoint)leftEye rightEye:(CGPoint)rightEye mouse:(CGPoint)mouse
{
    if (context) {
        CGContextClearRect(context, self.bounds) ;
    }
    context = UIGraphicsGetCurrentContext();
    CGContextAddRect(context, face);
    CGContextAddEllipseInRect(context, CGRectMake(leftEye.x - 1 , leftEye.y - 1 , 2 , 2));
    CGContextAddEllipseInRect(context, CGRectMake(rightEye.x - 1 , rightEye.y - 1 , 2 , 2));
    CGContextAddEllipseInRect(context, CGRectMake(mouse.x - 1 , mouse.y - 4 , 2 , 8));
    [[UIColor greenColor] set];
    CGContextSetLineWidth(context, 2);
    CGContextStrokePath(context);
}

@end
