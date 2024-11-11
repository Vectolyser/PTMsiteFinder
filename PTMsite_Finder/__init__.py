import pandas as pd
import re
import time

__name__ = 'PTMsite_Finder'


def extract_middle_columns(input_df, left_column_name='Precursor.Id', right_column_name='PTM.Site.Confidence'):
    # 读取 TSV 文件
    df = input_df

    # 获取列名列表
    column_names = df.columns.tolist()

    # 检查 A 列和 B 列是否在列名中
    if left_column_name in column_names and right_column_name in column_names:
        # 获取 A 列和 B 列的索引
        a_index = column_names.index(left_column_name)
        b_index = column_names.index(right_column_name)

        # 提取 A 列和 B 列之间的列名
        result_columns = column_names[a_index + 1:b_index]

        # 返回一个可被 pd.read_csv 读取的数组格式
        return pd.Series(result_columns)
    else:
        raise ValueError("One of the specified columns does not exist.")


def converting_rules(x, num):
    if num == 'STY(UniMod:21)':
        Mod = '(UniMod:21)'
    if num == 'K(UniMod:121)':
        Mod = '(UniMod:121)'
    if num == 'K(UniMod:1)':
        Mod = '(UniMod:1)'
    if num == 'N(UniMod:7)':
        Mod = '(UniMod:7)'
    if num == 'K(UniMod:92)':
        Mod = '(UniMod:92)'
    if num == 'STY(UniMod:21)':
        x = x.replace(f'S{Mod}', '1')
        x = x.replace(f'T{Mod}', '2')
        x = x.replace(f'Y{Mod}', '3')
        x = re.sub('\(UniMod:\d+\)', '', x)
    elif num == 'K(UniMod:121)':
        x = x.replace(f'K{Mod}', '4')
        x = re.sub('\(UniMod:\d+\)', '', x)
    elif num == 'K(UniMod:1)':
        x = x.replace(f'K{Mod}', '6')
        x = re.sub('\(UniMod:\d+\)', '', x)

    elif num == 'N(UniMod:7)':
        x = x.replace(f'N{Mod}', '5')
        x = re.sub('\(UniMod:\d+\)', '', x)
    elif num == 'K(UniMod:92)':
        x = x.replace(f'K{Mod}', '7')
        x = re.sub('\(UniMod:\d+\)', '', x)
    else:
        mod_str = num
        x = x.replace(f'{mod_str}', '8')
        x = re.sub('\(UniMod:\d+\)', '', x)
    return x


def timing_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # 记录开始时间  # 调用原函数
        func(*args, **kwargs)  # 调用函数
        end_time = time.time()  # 记录结束时间
        duration = end_time - start_time  # 计算执行时间
        return duration

    return wrapper


class UserStoppedException(Exception):
    """表示用户主动停止程序运行的异常。"""
    pass
