# 第五章：非约束性排序

本章内容简介：

- 5.1 Objectives
- 5.2 Ordination Overview
- 5.3 Principal Component Analysis (PCA)
- 5.4 Correspondence Analysis (CA)
- 5.5 Principal Coordinate Analysis (PCoA)
- 5.6 Nonmetric Multidimensional Scaling (NMDS)
- 5.7 Hand-Written PCA Ordination Function

## 5.1 Objectives

## 5.2 Ordination Overview

### 5.2.1 Multidimensional Space

### 5.2.2 Ordination in Reduced Space

## 5.3 Principal Component Analysis (PCA)

### 5.3.1 Overview

### 5.3.2 PCA of the Environmental Variables of the Doubs River Data Using rda()

#### 5.3.2.1 Preparation of the Data

#### 5.3.2.2 PCA of a Correlation Matrix

#### 5.3.2.3 Extracting, Interpreting and Plotting Results from a vegan Ordination Output Object

#### 5.3.2.4 Projecting Supplementary Variables into a PCA Biplot

#### 5.3.2.5 Projecting Supplementary Objects into a PCA Biplot

#### 5.3.2.6 Combining Clustering and Ordination Results

### 5.3.3 PCA on Transformed Species Data

#### 5.3.3.1 Application to the Hellinger-Transformed Fish Data

#### 5.3.3.2 Passive (post hoc) Explanation of Axes Using Environmental Variables

### 5.3.4 Domain of Application of PCA

### 5.3.5 PCA Using Function PCA.newr()

### 5.3.6 Imputation of Missing Values in PCA

## 5.4 Correspondence Analysis (CA)

### 5.4.1 Introduction

### 5.4.2 CA Using Function cca() of Package vegan

#### 5.4.2.1 Running the Analysis and Drawing the Biplots

#### 5.4.2.2 Projection of Sulementar Sites or Secies in a CA Bilot

#### 5.4.2.3 Post hoc Curve Fitting of Environmental Variables

#### 5.4.2.4 Reordering the Data Table on the Basis of an Ordination Axis

### 5.4.3 CA Using Function CA.newr()

### 5.4.4 Arch Effect and Detrended Correspondence Analysis (DCA)

### 5.4.5 Multiple Correspondence Analysis (MCA)

#### 5.4.5.1 MCA on the Environmental Variables of the Oribatid Mite Data Set

## 5.5 Principal Coordinate Analysis (PCoA)

### 5.5.1 Introduction

### 5.5.2 Application of PCoA to the Doubs Data Set Using cmdscale() and vegan

### 5.5.3 Application of PCoA to the Doubs Data Set Using pcoa()

## 5.6 Nonmetric Multidimensional Scaling (NMDS)

### 5.6.1 Introduction

### 5.6.2 Application to the Doubs Fish Data

### 5.6.3 PCoA or NMDS?

## 5.7 Hand-Written PCA Ordination Function

