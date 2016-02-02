//
//  ViewController.m
//  sample_face_track
//
//  Created by makun on 11/23/14.
//  Copyright (c) 2015 SenseTime.com . All rights reserved.
//

#import "ViewController.h"
#import <AVFoundation/AVFoundation.h>
#import "far.h"
#import "CanvasView.h"

@interface ViewController () <AVCaptureVideoDataOutputSampleBufferDelegate>

@property (nonatomic , strong) AVCaptureVideoPreviewLayer *captureVideoPreviewLayer ;

@property (nonatomic) cv_handle_t hTracker;

@property (nonatomic , strong) CanvasView *viewCanvas ;

@end

@implementation ViewController

- (void)viewDidLoad {
    
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    self.view.backgroundColor = [UIColor blackColor] ;
    
    self.hTracker = cv_face_create_tracker();
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (void)dealloc
{
    cv_face_destroy_tracker(self.hTracker);
}

- (void)viewDidAppear:(BOOL)animated
{
    [super viewDidAppear:animated];
    
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
    uint8_t *baseAddress = CVPixelBufferGetBaseAddress(pixelBuffer);
    
    cv_face_t *pFaceRectID = NULL ;
    int iCount = 0;
    
    int iWidth  = (int)CVPixelBufferGetWidth(pixelBuffer);
    int iHeight = (int)CVPixelBufferGetHeight(pixelBuffer);
    
    cv_result_t iRet = CV_OK;
    
    iRet = cv_face_track(self.hTracker, baseAddress, CV_PIX_FMT_BGRA8888, iWidth, iHeight, iWidth * 4, CV_FACE_LEFT, &pFaceRectID, &iCount);
    
    if ( iRet == CV_OK && iCount > 0 ) {
        
        NSMutableArray *arrPersons = [NSMutableArray array] ;
        
        for (int i = 0; i < iCount ; i ++) {
            
            cv_face_t rectIDMain = pFaceRectID[i] ;
            
            NSMutableArray *arrStrPoints = [NSMutableArray array] ;
            
            cv_rect_t rect = rectIDMain.rect ;
            
            CGRect rectFace = CGRectMake(rect.top , rect.left , rect.right - rect.left, rect.bottom - rect.top);
            
            NSMutableDictionary *dicPerson = [NSMutableDictionary dictionary] ;
            [dicPerson setObject:arrStrPoints forKey:POINTS_KEY];
            [dicPerson setObject:NSStringFromCGRect(rectFace) forKey:RECT_KEY];
            
            [arrPersons addObject:dicPerson] ;
        }
    
        dispatch_async(dispatch_get_main_queue(), ^{
            [self showFaceLandmarksAndFaceRectWithPersonsArray:arrPersons];
        } ) ;
        
    } else {
        dispatch_async(dispatch_get_main_queue(), ^{
            [self hideFace];
        } ) ;
    }
    cv_face_release_tracker_result(pFaceRectID, iCount);
    CVPixelBufferUnlockBaseAddress(pixelBuffer, 0);
}


- (void) showFaceLandmarksAndFaceRectWithPersonsArray:(NSMutableArray *)arrPersons
{
    if (self.viewCanvas.hidden) {
        self.viewCanvas.hidden = NO ;
    }
    self.viewCanvas.arrPersons = arrPersons ;
    [self.viewCanvas setNeedsDisplay] ;
}

- (void) hideFace {
    if (!self.viewCanvas.hidden) {
        self.viewCanvas.hidden = YES ;
    }
}

@end
