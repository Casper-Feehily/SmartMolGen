# **SmartMolGen - 智能分子结构生成器**

## **项目简介 | Project Description**

**SmartMolGen** 是一个基于 **RDKit** 的智能化学分子生成工具，用户只需输入碳、氢、氧原子数量，程序即可推测并生成可能的分子结构，包括：

- **烷烃 (Alkanes)**、**烯烃 (Alkenes)**、**炔烃 (Alkynes)**
- **芳香化合物 (Aromatic Compounds) - 苯、酚类等**
- **醇 (Alcohols)、酚 (Phenols)、醚 (Ethers)**
- **醛 (Aldehydes)、酮 (Ketones)、羧酸 (Carboxylic Acids)、酯 (Esters)**
- **自动识别并生成可能的环状化合物 (Cyclic Compounds)**

该工具可用于 **有机化学学习、化学结构预测、分子建模等**，为研究人员和学生提供帮助。

---

## **安装依赖 | Installation Requirements**

本项目基于 **Python (>=3.7)** 开发，并依赖 **RDKit** 进行化学结构的处理。

### **安装步骤 | Installation Steps**

1. **安装 Python**（建议使用 Anaconda 或 Miniconda）
2. **安装 RDKit**（可选方式）

    - **使用 Conda 安装（推荐）**：
      ```bash
      conda install -c conda-forge rdkit
      ```  
    - **使用 pip 安装（可能需要额外配置）**：
      ```bash
      pip install rdkit
      ```  

3. **安装其他依赖**（如未安装可手动执行）：
   ```bash
   pip install itertools

## 使用方法

- 运行 main。py，然后输入原子数量：
   ```bash
   python main.py
   ```
