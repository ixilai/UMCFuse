# UMCFuse: A Unified Multiple Complex Scenes Infrared and Visible Image Fusion Framework

![Workflow](https://github.com/ixilai/UMCFuse/blob/master/liucheng.jpg)

UMCFuse is a **unified framework** for **infrared and visible image fusion (IVIF)** in **complex scenes**, including haze, rain, snow, overexposure, fire, and noisy environments. It effectively integrates multi-modal information, preserves fine details, suppresses noise, and enhances scene understanding for downstream tasks such as **object detection**, **semantic segmentation**, and **salient object detection**.

---

## 🌟 Features

- **Unified Framework for Complex Scenes**  
  Handles multiple challenging environmental conditions within a single model, eliminating the need for scene-specific tuning.

- **Transmission Map Guided Decomposition**  
  Estimates **light transmission** for robust separation of contrast and structure layers, improving fusion quality in adverse conditions.

- **Adaptive High-Frequency Denoising**  
  Preserves fine details while suppressing noise in high-frequency components.

- **Multi-Scale Low-Frequency Extraction**  
  Effectively captures energy information in high-contrast regions for robust structural fusion.

---

## 🌟 Contrast and Structure Layer Decomposition

UMCFuse leverages **transmission maps** to decompose images into contrast and structure layers, separating degraded pixels from original content.  

- **Rainy scenes:** Raindrops scatter and refract light, creating blurred regions.  
- **Snow scenes:** Snowflakes scatter light, reducing contrast and capturing snow features.  
- **Overexposed areas:** Excessive light produces nearly uniform regions, highlighting saturated pixels.  

**Benefit:** Separating degradation reduces interference during feature extraction and improves the preservation of useful details.

### Figure 2: Transmission Map and Layer Decomposition
![fig2](https://github.com/ixilai/UMCFuse/blob/master/fig2.png)

### Figure 3: More Transmission Map Examples
![fig3](https://github.com/ixilai/UMCFuse/blob/master/fig3.jpg)

---

## 🖼 Qualitative Results

UMCFuse consistently preserves infrared luminance, visible textures, and scene details while suppressing noise across complex conditions.

### Normal, Noise, and Overexposure Scenes
![fig5](https://github.com/ixilai/UMCFuse/blob/master/fig4.png)  
*Fusion comparison on normal, noisy, and overexposed scenes.*

### Haze, Rain, and Blur Scenes
![fig6](https://github.com/ixilai/UMCFuse/blob/master/fig5.png)  
*Fusion comparison on haze, rain, and blur scenes.*

### Snow, Noise+Rain, and Fire Scenes
![fig7](https://github.com/ixilai/UMCFuse/blob/master/fig6.png)  
*Fusion performance in snow, rain+noise, and fire conditions.*

---

## 📊 Quantitative Comparison

We quantitatively evaluated UMCFuse against 11 state-of-the-art methods across seven datasets using metrics like $Q_G$, $EN$, $Q_{CV}$, $SSIM$, and more.

- **Overall Comparison:** Compared to the average of 11 competing methods, UMCFuse improves **13.64%** in $Q_{MI}$, $Q_{NCIE}$, $Q_G$, $Q_{abf}$ and **34.06%** in $Q_{CV}$, $VIF$, $EN$, $SSIM$.

### Table 1: Overall Performance
![Tab1](https://github.com/ixilai/UMCFuse/blob/master/Tab1.png)

### Table 2: Overall Performance
![Tab2](https://github.com/ixilai/UMCFuse/blob/master/Tab2.png)

---

## 🏥 Medical Image Fusion Experiments

UMCFuse demonstrates **cross-domain generalization** on medical images (50 SPECT-MRI pairs from Harvard Medical School). Compared with CDDFuse, SwinFusion, PRRGAN, and MDHU, UMCFuse consistently ranks **top two across five metrics**, highlighting robust feature extraction capability.

### Table 3: Medical Image Fusion Performance
![Tab3](https://github.com/ixilai/UMCFuse/blob/master/Tab3.png)

---

## 🧪 Downstream Tasks

UMCFuse enhances feature representation for tasks such as:  

- **Object Detection:** Higher detection rates and confidence.  
- **Semantic Segmentation:** Improved segmentation accuracy under adverse conditions.  
- **Salient Object Detection:** Accurate target localization with fewer artifacts.  
- **Depth Estimation:** Better depth maps by reducing interference from noise and haze.

---
## 🚀 Usage

## 🚀 Usage

The following demonstrates how to run **UMCFuse** on infrared and visible image pairs using the provided MATLAB demo.

### 1️⃣ Prerequisites

- MATLAB R2023a or other
- Image Processing Toolbox
- GPU recommended for faster computation

### 2️⃣ Directory Structure

Choice your source images

### 3️⃣ Run the Demo

Open `demo.m` in MATLAB (or any script)
