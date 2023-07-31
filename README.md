目录

- [1. 新版本(基于 rpy2)](#1-新版本基于-rpy2)
- [2. 开发过程所学](#2-开发过程所学)
  - [2.1. 背景知识 and 统计知识](#21-背景知识-and-统计知识)
  - [2.2. R 相关](#22-r-相关)
    - [2.2.1. 经验](#221-经验)
    - [2.2.2. API](#222-api)
  - [2.3. python list](#23-python-list)
  - [2.4. pandas and numpy 相关](#24-pandas-and-numpy-相关)
    - [2.4.1. 经验](#241-经验)
    - [2.4.2. API](#242-api)
  - [2.5. sklearn 相关](#25-sklearn-相关)
  - [2.6. scipy 相关](#26-scipy-相关)
  - [2.7. anndata 相关](#27-anndata-相关)
  - [2.8. pyvinecopulib 相关](#28-pyvinecopulib-相关)
  - [2.9. Statsmodels 相关](#29-statsmodels-相关)
    - [2.9.1. 经验](#291-经验)
    - [2.9.2. 补充](#292-补充)
  - [2.10. pygam 相关](#210-pygam-相关)
  - [2.11. python 相关](#211-python-相关)
    - [2.11.1. python 底层基础](#2111-python-底层基础)
    - [2.11.2. python 机制](#2112-python-机制)
  - [2.12. GitHub](#212-github)
  - [2.13. 其他](#213-其他)
- [3. 老版本(不用 rpy2 版本)的相关记录](#3-老版本不用-rpy2-版本的相关记录)
  - [3.1. 未来的下一步功能:](#31-未来的下一步功能)
  - [3.2. Questions](#32-questions)
  - [3.3. scDesign3 相关](#33-scdesign3-相关)
    - [3.3.1. construct data](#331-construct-data)
    - [3.3.2. fit marginal](#332-fit-marginal)

---

## 1. 新版本(基于 rpy2)

---

## 2. 开发过程所学

### 2.1. 背景知识 and 统计知识

[Tutorial](https://songdongyuan1994.github.io/scDesign3/docs/articles/scDesign3-introduction-vignette.html)

[GAM 建模、平滑系数选择、R 函数](https://www.bilibili.com/video/BV1h8411s7a1/?p=5&spm_id_from=333.880.my_history.page.click)

[mgcv_R_documentation](https://rdrr.io/cran/mgcv/api/)

[GAM 平滑参数估计方法 GCV 和 REML 的比较](https://github.com/DistanceDevelopment/dsm/wiki/Why-is-the-default-smoothing-method-%22REML%22-rather-than-%22GCV.Cp%22%3F)，REML 和 GCV 都试图达到同样的目标：使模型中的平滑部分恰到好处，不过分平滑。已经证明，当样本量趋向无穷大时，GCV 将选择最佳的平滑参数（即预测误差较低）。在实际应用中，可能存在一个样本量，使得这个渐近结果成立，但我们中很少有人拥有那么大规模的数据。在较小（有限）的样本量下，GCV 可能会出现多个最小值，导致优化困难，因此估计的平滑参数更加变化。此外，GCV 往往会低估平滑（即平滑程度比应该的要高），因为它对过度拟合的惩罚较弱（这在基于预测误差的准则中是有意义的，您希望模型能够很好地拟合数据）。REML（以及 ML、边际似然）对过度拟合的惩罚更重，因此具有更明显的极值，导致更少的优化问题和估计的平滑参数更加稳定。

[Statsmodels GAM documentation](https://www.statsmodels.org/stable/gam.html),[Statsmodels GLM documentation](https://www.statsmodels.org/stable/glm.html)

[pygam documentation](https://pygam.readthedocs.io/en/latest/index.html)

在广义可加模型（Generalized Additive Models，GAM）中，Edof（Effective Degrees of Freedom，有效自由度）是一种度量模型复杂度的指标。它表示模型在拟合数据时使用的有效自由度数量。在 GAM 中，每个平滑项（spline term）的自由度决定了该项对模型的贡献和灵活性。较高的自由度可以使平滑项更加灵活，能够更好地适应数据的变化，但也容易导致过拟合。Edof 考虑了平滑项的自由度以及其他相关因素，提供了一个更准确的模型复杂度度量。

[model selcetion](https://arxiv.org/pdf/1306.2427.pdf)

### 2.2. R 相关

#### 2.2.1. 经验

- 如何从 singlecellexperiment object 转换成 h5ad 存储后，用 scanpy 读取 [参考](https://www.bilibili.com/read/cv8675277/)使用`sceasy`包
- SummerizedExperiment 集成从 singlecellexperiment object 中提取 coldata/assay 等的方法，对应`anndata.obs`, `anndata.layer`
- `is.factor`等函数比起`typeof`更适合判断`data.frame`一列的数据类型
- 并行计算用函数 `parallel::mcmapply` `BiocParallel::bpmapply` `pbmcapply::pbmcmapply`都有两类参数，一类参数是直接传进去的，一类是封装在 MoreArgs 中传进去的，前者在每一个进程中是可迭代式地一个个取出来的，因此在每个进程中是不同的，后者是每一个进程中都需要用到的部分

```r
parafunc(fun_name,arg1=...,arg2=...,...,MoreArgs=list())
```

- `tryCatch`和 python 中的`try...except...`一致，`withCallingHandlers`，可以在特定的代码块中捕获和处理异常、警告和消息
- `list`和`data.frame`数据类型，`list`可以认为是一个字典（如果`list`创建的时候不写明键名字`list(0,1,2)`，会自动用 1 开始的数字编号作为键名字），访问某一个键的方式可以是`$name`也可以是`list[["name"]]`。`data.frame`可以看成一个特殊的`list`，列名是键名，每一个值都是一个`vector`对象，因此在访问列名时也可以使用`$name`方法，也因此调用`length()`函数时，`data.frame`输出的是列数，即键数目，和对`list`计算得到的结果是一样的（也和 python 对 dict 调用 len 输出键数是一样的）。当然`data.frame`还可以使用很多切片索引等方法是`list`不行的
- `lapply`和`sapply`的区别：都是对 list 的每个键下的元素用相同的 function 进行操作，区别在于`lapply`返回一个对应的 list，而`sapply`会尝试返回一个更简单的结果，如 vector 或者 matrix，若返回矩阵，其中**每列**对应于输入 list 中的一个元素的结果。（感觉操作的实际上是任意可迭代对象）
- `formula`对象可以用`update.formula`进行修改，在`formula`用`-`代表删去一些变量

#### 2.2.2. API

- `mclust:Mclust`使用高斯混合模型进行聚类 [参考](https://mclust-org.github.io/mclust/reference/Mclust.html)
- `interaction`函数计算 factor 变量的交互结果（乘积）
- `dplyr`操作 dataframe 的包 `dplyr::mutate(df,name = func(column))`: 在原来的 dataframe 中加一列名为`name`，内容为`func(column)`的列
- `table`: 计算每个元素的出现次数
- `names`: 返回给的 vector dataframe 等的元素名称
- singlecellexperiment 对象的存储和加载分别使用`saveRDS(sce,file ='xxx.rds')`和`loadRDS(sce,file ='xxx.rds')`实现
- `grepl`正则表达式
- `mgcv::gam(s(varible_name,k=10,bs='cr'),method=REML)`bs 代表回归样条的[基函数类型](https://www.rdocumentation.org/packages/mgcv/versions/1.9-0/topics/smooth.terms)（cr 为三次样条函数，三次样条函数指的就是具有连续性、且一阶和二阶导连续的三阶分段多项式），[k](https://www.rdocumentation.org/packages/mgcv/versions/1.9-0/topics/choose.k)，method 指定了平滑参数的估计方式
- `mgcv`中 s,ti,te 的区别和指定方式，简单而言，ti 仅代表交互效应，te 代表又有交互效应又有非交互效应，`te(x1,x2)`与`s(x1)+s(x2)+ti(x1,x2)`大致类似，但是估计的平滑惩罚数目不同(2 vs 4)。而 s(x1,x2)和 te(x1,x2)是类似的，区别在于 s(x1,x2)假设 x1,x2 有相同的 scale，即 isotropic 各向同性，te(x1,x2)认为是各变量的 scale 不同，各向异性(tensor product ti and te 处理的就是各向异性的情况)。[参考](https://stats.stackexchange.com/questions/519433/gam-and-multiple-continuous-continuous-interactions-tensor-smooths)
- `all.vars(formula)`返回一个 vector，第一位置是 y 的 var，第二位置是 x 的 var

### 2.3. python list

- list \* k，即可把 list 中的元素复制 k 遍，和 list1+list2 异曲同工
- `set(list)`得到 list 中所有的 unique 内容
- `map(func,iterable)`得到一个经过 func 操作的 iterable
- `set().issubset()`检查是不是一个类列表对象的子集，注意**空列表也会被认为是子集**
- 利用字典构造式进行字典合并

```python
def func(x):
    return {x: 'test'}

# 生成一个包含字典的列表
result_list = [func(i) for i in range(5)]

# 使用字典推导式将字典拼接成一个大字典
merged_dict = {k: v for d in result_list for k, v in d.items()}
```

### 2.4. pandas and numpy 相关

#### 2.4.1. 经验

- 从其他 dataframe 赋值的新 dataframe，要用.copy()得到一个备份，否则会出问题
- 如果用下述代码进行 dataframe 赋值，需要`series`和`dataframe`有相同的`index`，否则匹配不上的情况下`dataframe`中`index`匹配不上的部分直接赋值为`nan`，`series`的内容就丢了

```python
# not recommend, in this case, pd.concat may be preferred
df['name'] = pd.Series()
# recommend
df['name'] = np.array()
# or
df['name'] = list
# but not generator
```

- pandas 复制一行最快的办法：用这一行的索引重复调用该行即可

```python
# 得到index行重复repeat_times次的索引
repeat_index = np.repeat(np.array(index),repeat_times)

# 用该索引提取dataframe中的这一行
new_df = df.iloc[repeat_index]
```

- 遍历每一列的数据类型的方法(对于要把所有数值变量都认为是数值的方法，关键在于`np.issubdtype(dtype, np.number)`的应用)

```python
import pandas as pd
import numpy as np

# 创建示例DataFrame
data = {
    'col1': [1, 2, 3, 4],
    'col2': ['a', 'b', 'c', 'd'],
    'col3': [1.1, 2.2, 3.3, 4.4],
    'col4': [5, 6, 7, 8],
}
df = pd.DataFrame(data)

# 获取每列的数据类型
dtypes = df.dtypes

# 创建分类列表
category_cols = []
numeric_cols = []

# 遍历每列的数据类型
for col, dtype in dtypes.items():
    if dtype == 'category':
        category_cols.append(col)
    elif np.issubdtype(dtype, np.number):
        numeric_cols.append(col)

# 输出结果
print("Category columns:", category_cols)
print("Numeric columns:", numeric_cols)
```

- pandas 处理 NAN 数据(简单版)

```python
# 检查数据中是否有NAN最快的办法
test.isnull().values.any()
# inplace = True，则直接在原dataframe上进行修改
## 所有NA替换为100
test.fillna(100,inplace=True)
## 指定列的NA替换为某个给定的值
test.fillna({'lineage1':10,'lineage2':100},inplace=True)
## 删除指定列为NA的行
test.dropna(subset=['lineage2'],inplace=True)
## 删除所有feature中有NA的observation
test.dropna(inplace=True)
```

- Pandas 中的 isin() 方法来判断一个 Series 中的值是否在一个可迭代对象中。isin() 方法返回一个布尔型的 Series，表示每个元素是否在给定的可迭代对象中。

```python
# 示例数据
data = pd.Series([1, 2, 3, 4, 5])
values = [3, 5, 7]

# 判断每个元素是否在 values 中
result = data.isin(values)

# 0    False
# 1    False
# 2     True
# 3    False
# 4     True
# dtype: bool

# e.g 用来反选index
dat_use = dat_use.loc[~dat_use.index.isin(removed_cell)]
```

#### 2.4.2. API

- pandas 支持`dtype('category')`对应 R 中的`factor`类型 [参考](https://zhuanlan.zhihu.com/p/384727931)

```python
pd.dataframe['col'] = pd.dataframe['col'].astype('category')
```

- [`pd.groupby`操作](https://zhuanlan.zhihu.com/p/388789630)

```python
for name,group in df.groupby('name'):
    name    # category
    group   # corresponding dataframe
```

- `pd.value_counts`函数统计各个元素的出现次数
- `df.values`返回 dataframe 中的数值，为 np.array
- `df.describe()`返回所有的数值变量的描述性统计量，类似`summary` in R
- `np.random.seed(0)`设置`numpy`中的随机数种子

### 2.5. sklearn 相关

- `sklearn.mixture.GussianMixture`为高斯混合模型
- [KDE 核密度估计，非参方法估计分布函数](https://scikit-learn.org/stable/modules/density.html)
- sklearn 实现数据记忆的方式就是将一个 model 封装成一个类，类内不需要设置私有变量，可参考
- sklearn 训练数据如果是`np.array()`则维度必须为 n\*k，不能是只有单维度，也因此如果传 pandas 数据必须是`dataframe`不能是`series`

### 2.6. scipy 相关

- [稀疏矩阵信息](https://zhuanlan.zhihu.com/p/306320459)，可以使用`toarray()`转换成一般的 ndarray 的样式

### 2.7. anndata 相关

- anndata 数据结构存储都是 scipy 稀疏矩阵形式
- anndata.layer 存储不同类型的 assay（即对 anndata.X 的各种 transform）

```python
adata.layers['log_transformed'] = np.log1p(adata.X)
```

### 2.8. pyvinecopulib 相关

pyvinecopulib 和 rvinecopulib 是 C++库 Vinecopulib 在 python 和 R 上的分别实现

```python
import pyvinecopulib as pv
d = length(x.columns)
# calculate the empirical distribution of the raw RV
## u for a uniform distribution RV
u = pv.to_pseudo_obs(x)
# construct vine copula
cop = pv.Vinecop(data = u)
# simu new unifrom distribution RV
u_simu = cop.simulate(n=n_cell)
x_simu = np.asarray([np.quantile(x[:, i], u_sim[:, i]) for i in range(0, d)])
## x_simu shape n_variate * n_observation
x_simu = x_simu.T
```

输入数据 data 可以是 np.array 也可以是 pd.Dataframe，返回值为 np.array，其特性与 rvinecopulib 相似，只有一个变量时也能拟合与生成分布和采样数据。

注意其输入数据必须是 0-1 之间的数值，因为本质上 copula 函数是对一组 0-1 间的均匀分布数据进行计算，随机变量的 CDF 即为 0-1 均匀分布。但是 pyvinecopulib 不会像 rvinecopulib 一样给 raw data 自动做好转换为 CDF 的对应取值再做拟合，需要添加一步得到，并且要将拟合结果返回到原本的数据。

### 2.9. Statsmodels 相关

#### 2.9.1. 经验

- 使用`GLMGam.from_formula`构建模型会更好，线性的部分用 formula 的形式给出，其会将 dataframe 中 category 的部分直接转换为 dummy 编码并且自动添加截距项。如果直接用 y 和 x 对应的 dataframe 输入（用 endog 和 exog）则两者都需要手动操作，很不方便。`from_formula`由 package`patsy`提供支持([documentation](https://patsy.readthedocs.io/en/latest/))，`patsy`本身支持根据 mgcv 格式产生 spline 项，但是不支持惩罚项等的计算，因此无法投入应用。
- 超参数的选择可以指定，[`select_penweight`](https://www.statsmodels.org/stable/generated/statsmodels.gam.generalized_additive_model.GLMGam.select_penweight.html#statsmodels.gam.generalized_additive_model.GLMGam.select_penweight)，方法默认 aic，也可以选 gcv 等方法

```python
# 先大致指定建模方法
gam_bs = GLMGam.from_formula('YES ~ GDP + AGE + test',data = data.data, smoother=bs)
gam_bs.fit()
# 挑选超参数(平滑参数的选择)
alpha = gam_bs.select_penweight()[0]
# 指定完整的模型
model = GLMGam.from_formula('YES ~ GDP + AGE + test',data = data.data, smoother=bs,alpha=alpha)
# 模型参数拟合
result = model.fit()
```

- GLMGam 类使用 fit 方法后会返回 GLMGamResults 类，该类有两种获取 predict value 的方法，一种是`predict()`，直接返回预测值，为`np.array()`，另一种是`get_prediction().summary_frame()`前者返回一个 predcit result 类，后者用 dataframe 形式呈现结果，除了包含预测值外，还有对应的标准误差估计和置信区间

#### 2.9.2. 补充

[GAM 示例](https://www.statsmodels.org/stable/gam.html) [支持的 distribution](https://www.statsmodels.org/stable/glm.html) [model API](https://www.statsmodels.org/stable/generated/statsmodels.gam.generalized_additive_model.GLMGam.select_penweight.html#statsmodels.gam.generalized_additive_model.GLMGam.select_penweight)，`smoother`是样条函数选择，需要调用别的类来指定，在其中进一步规定样条函数自由度，`alpha`是平滑参数选择，可用`select_penweight()`得到。`k`代表了将一个拟合区间细分成 k+1 个小区间，每个区间使用一个 spline 去拟合，`k`与 spline 类型和要估计的参数`df`有换算关系。具体关系见[此](https://stats.stackexchange.com/questions/517375/splines-relationship-of-knots-degree-and-degrees-of-freedom#:~:text=In%20essence%2C%20splines%20are%20piecewise,degree%203%20and%20so%20on.)。

### 2.10. pygam 相关

- `model.gridsearch(X,y)`自动找到最好的超参数，[progress 是否显示进度条;用 lam 参数可以指定 lam 的搜索空间;objective 为搜索的优化目标](https://pygam.readthedocs.io/en/latest/api/gam.html)
- `model.sample`和`model.predict`的区别，sample 方法用于从模型的后验分布中生成样本。它可以用于生成模型的不确定性估计，以及对未观测数据的采样预测。sample 方法可以接受一个输入特征矩阵，并为每个样本生成一组预测。它返回一个 numpy 数组，其中每一行代表一个样本的预测。predict 方法用于对给定的输入进行预测。它接受一个输入特征矩阵，并为每个样本生成预测值。与 sample 方法不同，predict 方法不考虑模型的不确定性，它仅返回点估计的预测结果。它也返回一个 numpy 数组，其中每一行代表一个样本的预测值。
- AICc = AIC + 2p(p + 1)/(n − p − 1)

### 2.11. python 相关

#### 2.11.1. python 底层基础

- [如何组织 python 项目](https://zhuanlan.zhihu.com/p/335347908#:~:text=1%20%E7%BB%84%E7%BB%87%E7%BB%93%E6%9E%84%201%20Python%E9%A1%B9%E7%9B%AE%E7%9A%84%E7%BB%84%E7%BB%87%E7%BB%93%E6%9E%84%E4%B8%BB%E8%A6%81%E7%94%B1%20%E5%8C%85--%E6%A8%A1%E5%9D%97--%E7%B1%BB%20%E8%BF%99%E4%BA%9B%E9%83%A8%E5%88%86%E6%9E%84%E6%88%90%E3%80%82%202%20%E5%8C%85%EF%BC%9A,%E4%B8%80%E4%B8%AA%E6%A8%A1%E5%9D%97%E4%B8%8B%E9%9D%A2%E4%B9%9F%E5%8F%AF%E4%BB%A5%E5%8C%85%E5%90%AB%E5%A4%9A%E4%B8%AA%E7%B1%BB%E3%80%82%20%E6%A8%A1%E5%9D%97%E4%B8%8B%E9%9D%A2%E4%B9%9F%E5%8F%AF%E4%BB%A5%E7%9B%B4%E6%8E%A5%E5%86%99%E5%87%BD%E6%95%B0%E3%80%82%204%20%E7%B1%BB%EF%BC%9A%20%E7%B1%BB%E4%B8%AD%E5%AE%9A%E4%B9%89%E4%BA%86%E5%8F%98%E9%87%8F%E5%92%8C%E5%87%BD%E6%95%B0%E3%80%82%205%20%E5%87%BD%E6%95%B0%EF%BC%9A%20%E7%94%A8%E4%BA%8E%E5%AE%9E%E7%8E%B0%E7%89%B9%E5%AE%9A%E5%8A%9F%E8%83%BD%EF%BC%8C)
- `pip -V`查看现在使用的 pip 对应的下载位置和 python 解释器版本
- 线程与进程
  - [进程与线程的区别](https://www.bilibili.com/video/BV19b411L7y1/?spm_id_from=333.337.search-card.all.click&vd_source=44ff757ed2fdadeadf410d10bde17c2c)进程是分配资源的操作，是一个资源的综合，一个进程中是线程在进行操作，一个进程中的多线程，用一块资源来实现，实现多任务；而多进程则是多划了几块资源，每个进程只有一个主线程，每个线程都用自己的进程资源实现多任务
  - [python 的进程和线程的编程](https://www.bilibili.com/video/BV1jW411Y7pv?p=2&vd_source=44ff757ed2fdadeadf410d10bde17c2c)
  - [`multiprocessing`包大致功能介绍及编程](https://blog.csdn.net/sqzhao/article/details/120732881)
  - [如何得到子进程的返回值](https://blog.csdn.net/u013510614/article/details/109852590)

```python
# multiprocess demo
import multiprocessing as mp

# The function for the subprocess to call
def test(i):
    i += 1
    print(f'The {i}th process has been done')
    return i

if __name__ == '__main__':
    # Pool default for create the same number of subprocess as the cpu core
    pool = mp.Pool()
    storage = []
    # call the function name and the args are loaded as a tuple
    for i in range(15):
        storage.append(pool.apply_async(test, (i,)))
    # stop adding subprocess
    pool.close()
    # wait for every subprocess to finish
    pool.join()
    # get the result from the subprocess object
    result = [i.get() for i in storage]

    print(result)
```

- python 与 cpython 当使用 C++库提供的 Python 接口时，不需要额外安装 C++解释器。Python 本身已经提供了 CPython 解释器，它是使用 C 语言实现的，并且是 Python 的默认解释器。
- Anaconda 相关和库管理相关

  - `conda init powershell`使得 anaconda 在 windows powershell 中会显示虚拟环境名称，如果遇到遇到无法加载 profile.ps1，[解决方案](https://zhuanlan.zhihu.com/p/452273123)
  - 新 python 环境安装 jupyter kernel 的方式

  ```bash
  conda activate newenv
  conda install ipykernel
  python -m ipykernel install --user --name newenv --display-name 'any name is OK, just kernel name'
  ```

  - 依赖管理问题
    - `conda env export > environment.yaml`命令导出虚拟环境中的所有包依赖。用 conda 下载的包用`pip freeze > requirements.txt`时不会正确显示在 pypi 上的路径，只会显示一个本地路径，该命令主要是为了导出 pip 的相关依赖。
      - 个人感觉如果是想要包最后能直接发布到 Pypi 上，那就都用 Pypi 来装实现基础功能的包，这样方便直接得到直接可以用 Pypi 安装的 requirements，运行`pip install -r requirements.txt`即可实现迁移。
      - 如果为了真个 conda 环境的可迁移性，则直接用 ymal 文件即可，在创建新的虚拟环境时使用`conda env create -f environment.yaml`即可重新配置相同的环境。
    - 在 conda 环境中得到纯 pip 的 requirements 感觉也没什么特别好的办法，要么从 conda 的 yaml 里将 pip 相关的复制出来统一格式，要么从 pip freeze 中去掉本地安装路径的部分(目前在`modify_requirements.py`中实现，并验证过应该可以完整复现 pip 的所有安装)。<b>conda 和 pip 产生的文件都是 utf-16 编码的，处理解码时需要注意</b>
    - 如果全都是用 pip 装包，额外装了依赖，要卸除这些依赖，可以考虑`pip install python3-pip-autoremove`后运行`pip3-autoremove package_name`进行包及其依赖的卸载工作

- 通过 C++库的 Python 接口，可以在 Python 环境中直接调用 C++库的功能，而无需安装额外的 C++解释器。可以将 C++代码编译为动态链接库（如.so 文件或.dll 文件），然后在 Python 中导入这些库，并使用其中的函数和类。Python 解释器会与这些动态链接库进行交互，以实现对 C++代码的调用和执行。

#### 2.11.2. python 机制

- [import 函数详解](https://www.bilibili.com/video/BV1K24y1k7XA/?spm_id_from=333.788&vd_source=44ff757ed2fdadeadf410d10bde17c2c) 在 import 一个 module 的时候 python 会进行以下步骤
  - 创建模块命名空间：Python 会为模块创建一个独立的命名空间，该命名空间用于存储模块中定义的变量、函数和类。这样可以避免不同模块之间的命名冲突。
  - 执行模块代码：Python 会<b>按顺序执行</b>模块中的代码，将定义的变量、函数和类加载到模块的命名空间中。这样，其他代码可以通过模块名来访问这些定义。
  - 返回模块对象：在导入过程完成后，import module 语句会返回一个表示模块的对象，该对象可以在后续的代码中使用。可以使用该对象访问模块中定义的内容，例如变量、函数和类。
- [`if __name__ == '__main__'` 详解](https://blog.csdn.net/heqiang525/article/details/89879056)
- [函数指定传参的数据类型](https://blog.csdn.net/qq_42327755/article/details/87196150)
- [Python 类的私有公有变量设置](https://zhuanlan.zhihu.com/p/30223570)
- 判断 None Type 时，必须要用`variable is None`进行判断
- python warning 相关

```python
import warnings

# simplefilter 和 filterwarnings 可以认为前者是后者的简化版本
warnings.simplefilter('ignore')

def divide_numbers(a, b):
    if b == 1:
        warnings.warn(f"The {a}th divide is nonsense", RuntimeWarning)
        warnings.warn(f'Another warning for {a}',UserWarning)

## 两种warning都不会出现
divide_numbers(10, 1)

# 控制某个局部产生warning，并且抓住warning的信息(此时不会有warning的显示)
with warnings.catch_warnings(record =True) as warninglist:
    warnings.filterwarnings('default',category=RuntimeWarning)
    ## 只会出现runtime warning
    divide_numbers(9, 1)

# 重设warning filter
warnings.resetwarnings()
## 两种warning都会出现
divide_numbers(8, 1)

# 展示捕获的warning信息
for warning in warninglist:
    print(warning.message)
```

- 将脚本工作路径设置为当前脚本所在位置的方法

```python
import os
os.chdir(os.path.dirname(__file__))
```

- python 计时

```python
import time
start_time = time.time()
end_time = time.time()
running_time = start_time - end_time
```

- python 类型判断：`isinstance` or `type`，[建议使用前者](https://www.jianshu.com/p/7ef549503c93)
- python 在声明函数时如果可以接受不同的输入类型的话可以使用如下方式，[参考 1](https://www.bilibili.com/video/BV11Z4y1h79y/?spm_id_from=333.788.recommend_more_video.-1&vd_source=44ff757ed2fdadeadf410d10bde17c2c)，[参考 2](https://www.bilibili.com/video/BV16Y411K73d/?p=23&spm_id_from=pageDriver)，编译时通过标注类型实现编译时报错而非运行时报错

```python
from typing import Union
def toy(x:Union[list,str,int]) -> None:
    pass
```

- [python 装饰器](https://www.bilibili.com/video/BV1Gu411Q7JV/?spm_id_from=333.337.search-card.all.click&vd_source=44ff757ed2fdadeadf410d10bde17c2c)，本质上就是一个可以以函数为输入以函数为输出的函数，以`@decorator`形式写在装饰的函数的上方
- python 中可以通过`*args`和`**kwargs`进行**不定数量的参数**的传参，前者是以`tuple`的形式传入，后者以`dict`的形式传入

### 2.12. GitHub

[git 标准工作流程](https://www.bilibili.com/video/BV19e4y1q7JJ/?spm_id_from=333.337.search-card.all.click&vd_source=44ff757ed2fdadeadf410d10bde17c2c)

[获得 github token](https://blog.csdn.net/edsoki/article/details/126913354)

[git 操作介绍(菜鸟)](https://www.runoob.com/git/git-install-setup.html)

[git 版本管理与回退](https://www.bilibili.com/video/BV1ne4y1S7S9/?spm_id_from=333.337.search-card.all.click&vd_source=44ff757ed2fdadeadf410d10bde17c2c)，可以直接用 tag 名称（也可以用 harsh 值作为索引）作为回退，hard 会改变本地，soft 会只改变 commit，mixed 改变 commit 和 add。建议在版本回退之前先`git checkout -b`创建新分支在新分支上回退版本，避免丢失重要的东西。

[.gitignore 设置](https://cloud.tencent.com/developer/article/1640679)

```bash
# ignore all __pychche__
**/__pycache__
```

- 如果想要原本被跟踪的内容不再被跟踪，可以`git rm -r --cached .`把这些文件修改为不被跟踪的状态
- `.gitignore`和直接修改`.git/info/exclude`都能实现对文件的忽略，前者是对所有项目成员都有效，当所有项目成员都不需要版本控制一些文件的时候就用`.gitignore`，一般`.gitignore`也会上传到开源 repo 中，但是如果本地工作路径中有一些文件不需要做版本管理，那么直接按照`.gitignore`的格式将这些文件添加到`.git/info/exclude`中就好了，[参考 1](https://blog.csdn.net/qq_28505809/article/details/124946867)，[参考 2](https://www.jianshu.com/p/f42254bc3ffb)

[本地 github 传不上去的问题解决](https://blog.csdn.net/Xminyang/article/details/124837086)貌似是一个由于代理设置导致的问题

### 2.13. 其他

[VSCode 图标含义](https://www.zhihu.com/question/370258254)

[简单代码在线测试工具](https://app.datacamp.com/)

[VSCode 中使用 Black Formatter 和 isort](https://medium.com/mlearning-ai/python-auto-formatter-autopep8-vs-black-and-some-practical-tips-e71adb24aee1)，如果需要调整参数设置，可以在 arg 中设置，如`-l 100`设置每行最多的字符数；另外如果需要强制换行的，在最后多加一个`,`即可，formatter 会自动形成换行

---

## 3. 老版本(不用 rpy2 版本)的相关记录

### 3.1. 未来的下一步功能:

- [ ] 无法实现`fit marginal`中的`use bam`功能
- [ ] 边际分布函数族只包括 binomial，gaussian，poison，negative binomial，两种 zeroinflated 分布没做
- [ ] `fit marginal`使用的`statsmodels`的所有变量只能选择一种样条基函数(BSplines/CyclicCubicSplines)进行拟合
- [ ] 平滑超参数选择只能使用 gcv 方式选择
- [ ] 目前不支持`gamlss`模式
- [ ] 目前不支持`bam`模式
- [ ] 目前不支持 fit smooth surface，即无法考虑空间转录组数据(受限于`statsmodels`提供的功能有限)

### 3.2. Questions

- [x] vine copula 需要我手动去进行数据标准化吗(在 pyvinecopulib 中需要手动进行调整)
- [x] copulas package vine 在高于 3.8 的版本的应用
- [x] 用哪一种 criteria 进行超参数选择，默认是 AIC，需要换成 GCV 吗还有就是 CV 和 BIC (选 gcv 方法)

### 3.3. scDesign3 相关

#### 3.3.1. construct data

- [x] 一个 numeric 的 varible 就 fit kernel density estimation
- [x] pesudotime NA 是正常的，但是要 check 一下会不会报错(<b>会报错，需要解决，方案为直接退出程序让用户检查</b>)
- [x] warning 信息变为可选择性弹出
- [x] 尝试用 pyvinecopula 来 fit vine copula(<b>优先考虑使用 pyvinecopulib，该包与 R 中的 rvinecopulib 是相同的一个 C++库 vinecopulib 支持</b>)

#### 3.3.2. fit marginal

- mu sigma 都是只有~之后的 formula
- use bam bam 就是 gam，只是 for big data
- sigma 部分可以先留空，因为现在不用
- log 是记录 warning
- 183-210 去掉一部分不要的 cell，某一个 category 中如果有 gene 没有一个 cell 表达，会把这些 cell 去掉（而不是把 gene 去掉）
- mgcv 是 gam 的 package name
- `statsmodels`和`pygam`共性问题是支持的样条函数都比较少

|            `statsmodels`            |                 `pygam`                  |
| :---------------------------------: | :--------------------------------------: |
|           支持所有函数族            |                不支持 nb                 |
|       所有变量只能用一个样条        |        不同变量可以使用不同的样条        |
|              bs 和 cc               |     ps 和 cp(不知道是不是真的实现了)     |
| 不支持 te，不支持一切复杂一点的操作 | 支持 te，无法实现 mgcv 中 s(x1,x2)的效果 |
|             不支持 bam              |                不支持 bam                |

- 输出数据结构

```python
{
    'gene_name':{
        'fit': marginal_dis_model,
        'warning': {'warning_num': {
            'func_name':
            'type':
            'message':}}
        'time': fit_running_time,
        'removed_cell':
    }
}
```

- [x] 需要给用户调整样条函数的基函数和阶数的接口
- [x] 所有的 warning 返回给用户有一个 parameter
- [x] 记录时间，设置一个 para，如果用户需要就返回
- [x] patsy 无法 pickle，因此需要先更改数据阵的存储形式(用 patsy 生成数据阵后转换为 pandas 即可)再输入模型
- [ ] select penalty 运行速度奇慢无比，目前直接不做这个 select 了(调整一下优化方式，以后再选一个合适的算法，可能要大改包内的结构)
- [ ] statsmodels 需要广泛测试一下几个函数族能不能正确拟合出来，可能有一些会有问题，比如 binomial
- [ ] 加一个接口让用户可以调整一下各个 family 的参数设置

测试：

- [ ] alpha 和 gamma 的差距和区别（如何调出 R 中的 gamma 的数值？）
- [x] 比较 aic 的差距(自己和自己不带搜索，和 mgcv)
- [x] mgcv 中 negative binomial 的计算方法(在 statsmodels 中 $var=\mu+\alpha\mu^2$，在 mgcv 中 $var=\mu+\mu^2/\theta$，但是在 mgcv 中会优化这个 theta 参数使得尽量 dispersion 为 1，在 statsmodels 中是一个超参数)
- [x] mgcv 中的 select penal 的优化算法是什么（牛顿法）
- [x] 测试的 R 尽量与 python 对齐
- [x] 看看 statsmodel 超参数搜索为毛这么慢
