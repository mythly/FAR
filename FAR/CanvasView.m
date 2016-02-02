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
    [self drawPointWithPoints:self.arrPersons] ;
}

-(void)drawPointWithPoints:(NSArray *)arrPersons
{
    if (context) {
        CGContextClearRect(context, self.bounds) ;
    }
    context = UIGraphicsGetCurrentContext();
    
    for (NSDictionary *dicPerson in self.arrPersons) {
        if ([dicPerson objectForKey:POINTS_KEY]) {
            for (NSString *strPoints in [dicPerson objectForKey:POINTS_KEY]) {
                CGPoint p = CGPointFromString(strPoints) ;
                CGContextAddEllipseInRect(context, CGRectMake(p.x - 1 , p.y - 1 , 2 , 2));
            }
        }
        if ([dicPerson objectForKey:RECT_KEY]) {
            CGContextAddRect(context, CGRectFromString([dicPerson objectForKey:RECT_KEY])) ;
        }
    }

    [[UIColor greenColor] set];
    CGContextSetLineWidth(context, 2);
    CGContextStrokePath(context);
}

@end
