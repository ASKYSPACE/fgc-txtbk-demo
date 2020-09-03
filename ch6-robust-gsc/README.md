# 鲁棒变增益控制案例
本文件是本书第六章鲁棒变增益控制的配套程序。
其主要含有四个文件，其中:
- 主程序为 GainScheduledExample.m 
-  AirframeData.m  是飞行器与环境变量的参数设置
-  airframemodel.slx  是系统的模型
- airframecontrol.slx 是系统的控制回路。

本程序从‘GainScheduledExample.m‘进入，调用airframemodel.slx与airframecontrol.slx，其中:
- airframecontrol.slx会读取AirframeData.m的初始参数。
- 更加细致的案例可以参考mathwork中的’Tuning of Gain-Scheduled Three-Loop Autopilot‘案例，本文件是从其中获取的灵感与数据。

本文件中的‘airframemodel.slx’与‘airframecontrol.slx’是编写于matlab2019b版本中，用户在使用时需要使用相同或者更高版本的matlab，比如matlab2020.

