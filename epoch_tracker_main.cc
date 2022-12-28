/*
Copyright 2015 Google Inc. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <memory>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "core/file_resource.h"
#include "core/track.h"
#include "epoch_tracker/epoch_tracker.h"
#include "wave/wave.h"

// 这里是-help默认出来的信息。
const char* kHelp = "Usage: <bin> -i <input_file> "
    "[-f <f0_output> -p <pitchmarks_output> \\\n"
    "-c <correlations_output> "
    "-t "
    "-s "
    "-e <float> "
    "-x <float> "
    "-m <float> \\\n"
    "-u <float> "
    "-w <float> "
    "-a "
    "-d <debug_output_basename>] "
    "\n\n Help:\n"
    "-t enables a Hilbert transform that may reduce phase distortion\n"
    "-s suppress applying high pass filter at 80Hz "
    "(rumble-removal highpass filter)\n"
    "-e specifies the output frame interval for F0\n"
    "-x maximum f0 to look for\n"
    "-m minimum f0 to look for\n"
    "-u regular inter-mark interval to use in UV pitchmark regions\n"
    "-w set the cost for unvoiced segments\n"
    "   (def. 0.9, the higher the value the more f0 estimates in noise)\n"
    "-a saves F0 and PM output in ascii mode\n"
    "-d write diagnostic output to this file pattern\n"
"\nOutputs:\n"\
"The output files specified with -f and -p are Edinburgh Speech Tool\n"
"(EST) style files.  These are output as binary by default, ASCII with\n"
"the -a option.  The F0 values in the f0_output file are resampled with\n"
"the frame interval specified by option -e (default .005 s).\n"
"The unvoiced regions of the pitchmark file are filled with marks spaced\n"
"by the interval specified with -u (default .01 s).\n"
"\nIt is strongly recommended that the default high-pass filter be\n"
"applied to all input signals to remove DC and low-frequency\n"
"components.  For signals that have been recorded using close-talking\n"
"microphones, or those that have been subjected to various other\n"
"non-linear phase distortions in the studio, or in post-production, it\n"
"is often helpful to apply a Hilbert transform (-t option).  The .resid\n"
"file output when -d is specified can be examined to determine if the\n"
"voice pulses look anything like the classical glottal-flow derivative,\n"
"and the Hilbert transform enabled, or not.\n"
"\n"
"In the discussion below, the following notation is used:\n"
"Fs: sample rate\n"
"Fs0: sample rate of input file\n"
"nd: n-dimensional vector signal\n"
"ts: time-stamped vector signal; 1st ele. is the time of the sample in sec.\n"
"float, int16, etc.: atomic type of the data in the file\n"
"\nIf -d <name>is specified, the following raw binary output files are produced,\n"
"with the base path as specified in <name>, and the indicated extensions:\n"
".bestcorr 2d float ts the highest NCC value found at each residual peak\n"
".bprms 1d float Fs=500 RMS of 100-1000 Hz bandpassed input signal\n"
".f0ap 3d float ts F0 and NCC value found for each period\n"
".resid 1d float Fs=Fs0 LPC residual of conditioned input signal\n"
".nresid 1d float Fs=Fs0 LPC residual with local gain normalization\n"
".offsetp 1d float Fs=500 pseudo-probability that voicing is terminating\n"
".onsetp 1d float Fs=500 pseudo-probability that voicing is starting\n"
".pvoiced 1d float Fs=500 pseudo-probability that voicing is occurring\n"
".pcm 1d float Fs=Fs0 conditioned input signal\n"
".pmlab ASCII 'xlabel' format file of epoch marks\n"
".pvals 1d float Fs=Fs0 graded residual peak candidates\n";

Track* MakeEpochOutput(EpochTracker &et, float unvoiced_pm_interval) {
  std::vector<float> times;
  std::vector<int16_t> voicing;
  et.GetFilledEpochs(unvoiced_pm_interval, &times, &voicing);
  // https://cplusplus.com/doc/tutorial/classes/#pointers_to_classes
  // Similarly as with plain data structures, the members of an object can be accessed directly from a pointer by using the arrow operator (->).

  // Track* pm_track is a pointer to an object of class Track.
//
//    栈： 在函数内部声明的所有变量都将占用栈内存
//    堆： 这是程序中未使用的内存，在程序运行时可用于动态分配内存
  Track* pm_track = new Track; // creates a track of type track
//    https://www.twle.cn/l/yufei/cplusplus/cplusplus-basic-dynamic-memory.html
//    在 C++ 中，我们可以使用 new 运算符为给定类型的变量在运行时分配堆内的内存，然后返回所分配的空间地址, new data-type;
//
//    一旦不需要动态分配的内存，可以使用 delete 运算符来删除由 new 运算符分配的内存
// 比如我们需要定义一个指向 Track 类型的指针pm_track，然后请求内存，该内存在执行时被分配。


  pm_track->resize(times.size());
  for (int32_t i = 0; i < times.size(); ++i) {
    pm_track->t(i) = times[i];
    pm_track->set_v(i, voicing[i]);
  }
  return pm_track;
}


// declare cor as a pointer to a pointer to an int
Track* MakeF0Output(EpochTracker &et, float resample_interval, Track** cor) {
  std::vector<float> f0;
  std::vector<float> corr;
  if (!et.ResampleAndReturnResults(resample_interval, &f0, &corr)) {
    return NULL;
  }

  Track* f0_track = new Track;
  Track* cor_track = new Track;
  f0_track->resize(f0.size());
  cor_track->resize(corr.size());
  for (int32_t i = 0; i < f0.size(); ++i) {
    float t = resample_interval * i;
    // It is used with a pointer variable pointing to a structure or union.
    f0_track->t(i) = t;
    // The Dot(.) operator is used to normally access members of a structure or union.
    cor_track->t(i) = t;
    // The Arrow(->) operator exists to access the members of the structure or the unions using pointers.
    f0_track->set_v(i, (f0[i] > 0.0) ? true : false);
    cor_track->set_v(i, (f0[i] > 0.0) ? true : false);
    f0_track->a(i) = (f0[i] > 0.0) ? f0[i] : -1.0;
    cor_track->a(i) = corr[i];
  }
  *cor = cor_track;
  return f0_track;
}

// 这里会返回布尔值。
bool ComputeEpochsAndF0(EpochTracker &et, float unvoiced_pulse_interval,
			float external_frame_interval,
			Track** pm, Track** f0, Track** corr) {
  if (!et.ComputeFeatures()) {
    return false;
  }
  bool tr_result = et.TrackEpochs();
  et.WriteDiagnostics("");  // Try to save them here, even after tracking failure.
  if (!tr_result) {
    fprintf(stderr, "Problems in TrackEpochs");
    return false;
  }

  // create pm and f0 objects, these need to be freed in calling client.
  *pm = MakeEpochOutput(et, unvoiced_pulse_interval);
  *f0 = MakeF0Output(et, external_frame_interval, corr);
  return true;
}



//To compute F0 (pitch) and pitchmark (GCI) tracks and write them out as ASCII files:
//reaper -i /tmp/bla.wav -f /tmp/bla.f0 -p /tmp/bla.pm -a

int main(int argc, char* argv[]) {

//static const float kExternalFrameInterval = 0.005;
//static const float kInternalFrameInterval = 0.002;
//static const float kMinF0Search = 40.0;
//static const float kMaxF0Search = 500.0;
//static const float kUnvoicedPulseInterval = 0.01;
//static const float kUnvoicedCost = 0.9;
// 默认开启高通滤波，-s的话就会关闭
//static const bool kDoHighpass = true;
//  true 就是布尔值1
//static const bool kDoHilbertTransform = false;
//static const char kDebugName[] = "";


  int opt = 0;
  std::string filename;
  // 这里的f0_output是输出的f0的文件名字
  std::string f0_output;
  // 这里是GCI track文件
  std::string pm_output;

  std::string corr_output;
  bool do_hilbert_transform = kDoHilbertTransform;
  bool do_high_pass = kDoHighpass;
  float external_frame_interval = kExternalFrameInterval;
  float max_f0 = kMaxF0Search;
  float min_f0 = kMinF0Search;
  float inter_pulse = kUnvoicedPulseInterval;
  float unvoiced_cost = kUnvoicedCost;
  bool ascii = false;
  std::string debug_output;
  if (argc < 3) {
    fprintf(stdout, "\n%s\n", kHelp);
    return 1;
  }
  while ((opt = getopt(argc, argv, "i:f:p:c:htse:x:m:u:w:ad:")) != -1) {
    switch(opt) {
      case 'i':
        filename = optarg;
        break;
      case 'f':
        f0_output = optarg;
        break;
      case 'p':
        pm_output = optarg;
        break;
      case 'c':
        corr_output = optarg;
        break;
      case 't':
        do_hilbert_transform = true;
        break;
      case 's':
          // 这里决定了我们关不关高通滤波器
        do_high_pass = false;
        break;
      case 'e':
        external_frame_interval = atof(optarg);
        break;
      case 'x':
        max_f0 = atof(optarg);
        break;
      case 'm':
        min_f0 = atof(optarg);
        break;
      case 'u':
        inter_pulse = atof(optarg);
        break;
      case 'w':
        unvoiced_cost = atof(optarg);
        break;
      case 'a':
        ascii = true;
        break;
      case 'd':
        debug_output = optarg;
        break;
      case 'h':
        fprintf(stdout, "\n%s\n", kHelp);
        return 0;
    }
  }

  // Load input.
  Wave wav;
  if (!wav.Load(filename)) {
    fprintf(stderr, "Failed to load waveform '%s'\n", filename.c_str());
    return 1;
  }

  EpochTracker et;
  // void set_unvoiced_cost(float v) { unvoiced_cost_ = v; }
  et.set_unvoiced_cost(unvoiced_cost);
  // unvoiced_cost在这里是kUnvoicedCos等于0.9
  // 16 bit 的采样率，这里转成了16bit的数据，放在指针wave_datap里面
  int16_t* wave_datap = const_cast<int16_t *>(wav.data()->data());
  int32_t n_samples = wav.num_samples();
  float sample_rate = wav.sample_rate();
    // Prepare the instance for use.  Some sanity check is made on the
    // parameters, and the instance is reset so it can be reused
    // multiple times by simply calling Init() for each new input
    // signal. frame_interval determines the framing for some of the
    // internal feature computations, and for the periodic resampling of
    // F0 that will occur during final tracking result output.  min_
    // and max_f0_search are the bounding values, in Hz, for the F0
    // search.
    // NOTE: This Init method is DEPRECATED, and is only retained to
    // support legacy code.  IT MAY GO AWAY SOON.  This is NOT to be
    // used with ComputeEpochs().

    // 初始化对象et
    // The way to understand the “!” operator is to think that as a “not”,
    // so when you ask if(!programming) is a way to express “if it does not program do such”….
    // You could read it as "when there is no Init".
  if (!et.Init(wave_datap, n_samples, sample_rate,
	       min_f0, max_f0, do_high_pass, do_hilbert_transform)) {
    return 1;
  }
  if (!debug_output.empty()) {
    et.set_debug_name(debug_output);
  }
  // Compute f0 and pitchmarks.
  Track *f0 = NULL;  // // 初始化为 null 的指针f0
  Track *pm = NULL;
  Track *corr = NULL;
  if (!ComputeEpochsAndF0(et, inter_pulse, external_frame_interval, &pm, &f0, &corr)) {
    fprintf(stderr, "Failed to compute epochs\n");
    return 1;
  }

  // Save outputs.
  //   bool Save(const std::string &filename, bool ascii) const;
  // (pointer_name)->(variable_name)
  // https://www.geeksforgeeks.org/arrow-operator-in-c-c-with-examples/
  // Assigning value to Save(f0_output, ascii))


  // save 函数Save(const std::string &filename, bool ascii)
  if (!f0_output.empty() && !f0->Save(f0_output, ascii)) {
    delete f0;
    fprintf(stderr, "Failed to save f0 to '%s'\n", f0_output.c_str());
    return 1;
  }
  if (!pm_output.empty() && !pm->Save(pm_output, ascii)) {
    delete pm;
    fprintf(stderr, "Failed to save pitchmarks to '%s'\n", pm_output.c_str());
    return 1;
  }
  if (!corr_output.empty() && !corr->Save(corr_output, ascii)) {
    delete corr;
    fprintf(stderr, "Failed to save correlations to '%s'\n", corr_output.c_str());
    return 1;
  }
  delete f0;  // 释放 f0 所指向的内存
  delete pm;
  delete corr;
  return 0;
}
