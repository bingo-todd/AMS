# AMS
amplitude modulation spectrum
参考[AMS matlab工具包](http://ecs.utdallas.edu/loizou/AMS_Binary_Mask_Demos/demos.html)完成，舍弃了Matlab工具包中与AMS无关的部分。
AMS即幅度调制谱，计算过程如下：
-  分带滤波，使用butterworth滤波器（多个3阶带同滤波器+1个6阶高通）；
-  在每个频带内，对信号进行半波整流（也可以直接对信号取square）。之后是提取每个频带内的时域包络，这里跟平时的做法不同，**作者用matlab中的decimate函数将信号将采样至4000 Hz，将降采样之后的信号看作包络。**
-  在每个频带内，计算包络信号的STFT，帧长为128点，对应32ms（采样率为4000Hz），帧移为帧长的一半。傅立叶变换的点数为256，等效在每帧包络信号之后补128个0
-  将128个有效FFT系数压缩至15个。在0～400Hz范围内，设计15个均匀分布的三角窗，每个三角窗与FFT相乘即可得到15个系数。
Matlab工具包中将各频带的AM特征拼接成单维向量，该过程并不是将各频带的AM特征直接首尾拼接，而是将所有频带同一维度的特征拼在一起，

在Matlab中
ams=[feature_len,freq_channel_num]
```matlab
ams = ams';
reshape(ams,[],1);
```
