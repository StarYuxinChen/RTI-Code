This is a series of Matlab codes for my own reference. The codes is written by Stefan Nixon. 

mainSID.m是主代码，把原始的 LIF 灰度图像，变成一个“已经校正了激光、不均匀照明、条纹噪声和光学几何畸变的、可以当作相对浓度场来用的标量场”。
整个流程分成三块儿：
1. 标定与坐标变换：把相机坐标映射到一个“光线竖直、光路平行”的坐标系，方便在物理意义上沿激光传播方向做运算。
2. 背景与激光片的确定：用 4 张背景图（fIn2代表无染料有灯光；fIn3代表无染料无灯光；fIn4代表有染料有灯光；fIn5代表有染料无灯光）来估计：相机自身噪声；激光片空间分布；染料均匀时的响应（用于标度）。
3. 对每一帧实验 LIF 图像做：校正照明不均匀；选择性去条纹；考虑激光沿程衰减和染料的吸收；归一化到一个“0–1”的相对浓度。（使用correct_Plif_MI_new这个函数，包含使用小波变换和FFT）
Note：想要了解更多有关小波转换的信息可以看munch et al. 2009: "Stripe and ring artifact removal with combined wavelet — Fourier filtering"

其余函数有各自的功能：
1. dfi2mat：把 DigiFlow的图片(.dfi)文件变成 MATLAB 能算的图像，更具体地说转换成了一个Matlab结构体Struct。可以随时调取图像数据、时间、空间标度、标定信息等。
2. mapTo：用标定结果把图像转换到目标坐标系，含义：根据 x_map_file 和 y_map_file，把原始图像重采样到“光线平行坐标系”。
3. RemoveStripesVertical：竖直条纹去除滤波，含义：用“小波分解 + 频域滤波”的方法，把图像中的竖直条纹 artefacts 去掉。小波可以把图像分解到不同的空间尺度（大结构 vs 小结构）；在每个尺度里再做 FFT，可以识别“特定方向的周期性条纹”。
4. df_dfm_info:给定一个DigiFlow 生成的 .dfm 电影文件，读取文件头和帧目录，返回一个结构体 p，里面包含：这段视频的基本信息（谁做的、什么时候做的、多大、多长）、每一帧在文件中的位置、以及图像的大小等。
5. 这个 df_dfm_read 可以理解成：把 .dfm 电影文件里的若干帧，真正读成 MATLAB 里的图像数组的函数。
6. dfiRegion：按照 DigiFlow 里选的 extract_region，把 dfi 图像在 MATLAB 里裁出来，而且把“底部为原点”的坐标系换成 MATLAB 的“左上角为原点”。因为digiFlow和matlab的图片坐标系不一致。
7. dfiPlot：把 dfi2mat 读进来的 dfi 结构体，自动画成图：先画标量场（image），再叠加矢量场（vectors）
8. fillZeros：把矩阵里为 0 的点，用它上下左右非零邻居的平均值来“填补”，相当于做一个很简单的插值 / 填洞操作。
9. FiltWindow：对整张图像，用一个 window_size × window_size 的小框在上面滑过去，每次在这个框里做一点处理，然后把结果累加起来，最后平均一下，得到一张被“抹平/降噪”后的图片。
10. text2mat：把一个存成 .txt 的数字表格（比如坐标映射数组）读进来，按指定的行列数还原成 MATLAB 里的矩阵。

做实验需要拍摄的：
1. 4个背景图：fIn2代表无染料有灯光；fIn3代表无染料无灯光；fIn4代表有染料有灯光；fIn5代表有染料无灯光
2. 不对齐的点格：用于做
3. Streaks：用于生成 x_map_file 和 y_map_file，把原始图像重采样到“光线平行坐标系”。
Note：生成坐标系的代码在digiflow里
