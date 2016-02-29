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
    [self drawFace:CGRectFromString(self.strFace) leftEye:CGPointFromString(self.strLeftEye) rightEye:CGPointFromString(self.strRightEye) mouse:CGPointFromString(self.strMouse) info:self.strInfo] ;
}

-(void)drawFace:(CGRect)face leftEye:(CGPoint)leftEye rightEye:(CGPoint)rightEye mouse:(CGPoint)mouse info:(NSString*)info
{
    if (context) {
        CGContextClearRect(context, self.bounds) ;
    }
    context = UIGraphicsGetCurrentContext();
    
    
    NSMutableParagraphStyle* textStyle = NSMutableParagraphStyle.defaultParagraphStyle.mutableCopy;
    textStyle.alignment = NSTextAlignmentLeft;
    NSDictionary* textFontAttributes = @{NSFontAttributeName: [UIFont fontWithName: @"Helvetica" size: 16], NSForegroundColorAttributeName: UIColor.redColor, NSParagraphStyleAttributeName: textStyle};
    CGContextSaveGState(context);
    CGContextTranslateCTM(context, 60, self.frame.size.height - 10);
    CGContextRotateCTM(context, -(M_PI/2));
    [info drawInRect: CGRectMake(0, 0, 120, 120) withAttributes: textFontAttributes];
    CGContextRestoreGState(context);
    
    CGContextAddRect(context, face);
    CGContextAddEllipseInRect(context, CGRectMake(leftEye.x - 1 , leftEye.y - 1 , 2 , 2));
    CGContextAddEllipseInRect(context, CGRectMake(rightEye.x - 1 , rightEye.y - 1 , 2 , 2));
    CGContextAddEllipseInRect(context, CGRectMake(mouse.x - 1 , mouse.y - 4 , 2 , 8));
    
    [[UIColor greenColor] set];
    CGContextSetLineWidth(context, 2);
    CGContextStrokePath(context);
}

@end
