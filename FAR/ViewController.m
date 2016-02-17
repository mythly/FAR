//
//  ViewController.m
//  sample_face_track
//
//  Created by makun on 11/23/14.
//  Copyright (c) 2015 SenseTime.com . All rights reserved.
//

#import "ViewController.h"
#import <AVFoundation/AVFoundation.h>
#import <ImageIO/ImageIO.h>
#import "CanvasView.h"
#import "cv_face.h"

@interface ViewController () <AVCaptureVideoDataOutputSampleBufferDelegate>

@property (nonatomic , strong) CIContext *context;

@property (nonatomic , strong) AVCaptureVideoPreviewLayer *captureVideoPreviewLayer ;

@property (nonatomic , strong) CanvasView *viewCanvas ;

@end

@implementation ViewController

- (void)viewDidLoad {
    
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    self.view.backgroundColor = [UIColor blackColor] ;
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (void)dealloc
{
    cv_release();
}

- (void)viewDidAppear:(BOOL)animated
{
    [super viewDidAppear:animated];
    
    self.context = [CIContext contextWithOptions:nil];
    
    AVCaptureSession *session = [[AVCaptureSession alloc] init];
    
    session.sessionPreset = AVCaptureSessionPreset640x480;
    
    self.captureVideoPreviewLayer = [[AVCaptureVideoPreviewLayer alloc] initWithSession:session];
    self.captureVideoPreviewLayer.frame = CGRectMake( 0, 0, 480, 640 ) ;
    self.captureVideoPreviewLayer.position = self.view.center ;
    [self.captureVideoPreviewLayer setVideoGravity:AVLayerVideoGravityResizeAspectFill];
    [self.view.layer addSublayer:self.captureVideoPreviewLayer];
    
    self.viewCanvas = [[CanvasView alloc] initWithFrame:self.captureVideoPreviewLayer.frame] ;
    [self.view addSubview:self.viewCanvas] ;
    self.viewCanvas.backgroundColor = [UIColor clearColor] ;
    
    AVCaptureDevice *deviceFront ;
    
    NSArray *devices = [AVCaptureDevice devices];
    for (AVCaptureDevice *device in devices) {
        if ([device hasMediaType:AVMediaTypeVideo]) {
            
            if ([device position] == AVCaptureDevicePositionFront) {
                deviceFront = device;
            }
        }
    }
    
    NSError *error = nil;
    AVCaptureDeviceInput *input = [AVCaptureDeviceInput deviceInputWithDevice:deviceFront error:&error];
    if (!input) {
        // Handle the error appropriately.
        NSLog(@"ERROR: trying to open camera: %@", error);
    }
    AVCaptureVideoDataOutput * dataOutput = [[AVCaptureVideoDataOutput alloc] init];
    [dataOutput setAlwaysDiscardsLateVideoFrames:YES];
    [dataOutput setVideoSettings:[NSDictionary dictionaryWithObject:[NSNumber numberWithInt:kCVPixelFormatType_32BGRA] forKey:(id)kCVPixelBufferPixelFormatTypeKey]];
    
    dispatch_queue_t queue = dispatch_queue_create("bufferQueue", NULL);
    [dataOutput setSampleBufferDelegate:self queue:queue];
    
    [session beginConfiguration];
    if ([session canAddInput:input]) {
        [session addInput:input];
    }
    if ([session canAddOutput:dataOutput]) {
        [session addOutput:dataOutput];
    }
    [session commitConfiguration];
    
    [session startRunning];
}

-(void)captureOutput:(AVCaptureOutput *)captureOutput didOutputSampleBuffer:(CMSampleBufferRef)sampleBuffer fromConnection:(AVCaptureConnection *)connection {
    
    CVPixelBufferRef pixelBuffer = (CVPixelBufferRef)CMSampleBufferGetImageBuffer(sampleBuffer);
    CVPixelBufferLockBaseAddress(pixelBuffer, 0);
    uint8_t* baseAddress = CVPixelBufferGetBaseAddress(pixelBuffer);
    
    int iWidth  = (int)CVPixelBufferGetWidth(pixelBuffer);
    int iHeight = (int)CVPixelBufferGetHeight(pixelBuffer);
    unsigned char *gray = malloc(iWidth * iHeight);
    for (int i = 0; i < iHeight * iWidth; ++i)
        gray[i] = 0.299f * baseAddress[i * 4 + 2] + 0.587f * baseAddress[i * 4 + 1] + 0.114f * baseAddress[i * 4];
     
    cv_rect_t rect;
    float l = (iWidth < iHeight ? iWidth : iHeight) * 0.25f;
    rect.x = iWidth * 0.5f - l * 0.5f;
    rect.y = iHeight * 0.5f - l * 0.5f;
    rect.width = rect.height = l;
    
    NSLog(@"new frame %dx%d", iWidth, iHeight);
    if (cv_check()) {
        rect = cv_face_track(gray);
        NSLog(@"track at [(%.0f,%.0f) %.0fx%.0f]", rect.x, rect.y, rect.width, rect.height);
    }else {
        NSDictionary *opts = @{ CIDetectorAccuracy : CIDetectorAccuracyLow };
        CIDetector *detector = [CIDetector detectorOfType:CIDetectorTypeFace context:self.context options:opts];
        CIImage *ciImage = [CIImage imageWithCVPixelBuffer:pixelBuffer];
        opts = @{ CIDetectorImageOrientation : [NSNumber numberWithInt:6]};
        NSArray *features = [detector featuresInImage:ciImage options:opts];
        if ([features count] > 0) {
            CIFeature *f = [features objectAtIndex:0];
            rect.x = f.bounds.origin.x;
            rect.y = f.bounds.origin.y;
            rect.width = f.bounds.size.width;
            rect.height = f.bounds.size.height;
            cv_init(gray, iWidth, iHeight, rect);
            NSLog(@"detect at [(%.0f,%.0f) %.0fx%.0f]", rect.x, rect.y, rect.width, rect.height);
        }else
            NSLog(@"detect nothing\n");
    }

    CGRect rectFace = CGRectMake(iHeight - rect.height - rect.y, rect.x, rect.height, rect.width);
    dispatch_async(dispatch_get_main_queue(), ^{
        [self showFace:rectFace];
    } ) ;
    
    free(gray);
    CVPixelBufferUnlockBaseAddress(pixelBuffer, 0);
}


- (void) showFace:(CGRect)rectFace
{
    self.viewCanvas.strFace = NSStringFromCGRect(rectFace);
    [self.viewCanvas setNeedsDisplay] ;
}

@end
