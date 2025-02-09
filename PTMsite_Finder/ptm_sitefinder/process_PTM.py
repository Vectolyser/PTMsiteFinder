from . import extract_middle_columns


def ptm_sum(input_df, output_file,
            a_column_name='Precursor.Id',
            b_column_name='PTM.Site.Confidence'):
    # 读取 TSV 文件
    df = input_df

    df = df[df['Modified.Amino.Acid'].notnull() & (df['Modified.Amino.Acid'] != '')]

    # 将 Protein.Group_Position 列切分，并展开成多行
    df['Protein.Group_Position'] = df['Protein.Group_Position'].str.split(',')
    df['Sequence.Window'] = df['Sequence.Window'].str.split(',')
    expanded_df = df.explode(['Protein.Group_Position', 'Sequence.Window'])

    # 创建 PTM 列
    expanded_df['Protein.PTMsite'] = expanded_df['Protein.Group_Position']

    # 提取中间列
    extracted_columns = extract_middle_columns(df, a_column_name, b_column_name)

    # 进行分组，并对提取的列进行求和
    grouped_df = expanded_df.groupby('Protein.PTMsite').agg(
        {**{col: 'first' for col in df.columns if col not in ['Protein.PTMsite']},  # 保留第一条记录的其他列
         **{col: 'sum' for col in extracted_columns}, 'Sequence.Window': 'first'}  # 对提取列进行求和
    ).reset_index()

    pick_columns = \
        ['Protein.PTMsite', 'Protein.Names', 'Genes', 'Sequence.Window', 'PTM.Site.Confidence'] + \
        extracted_columns.tolist()

    grouped_df = grouped_df[pick_columns]

    # 将结果保存为新的 TSV 文件
    grouped_df.to_csv(output_file, sep='\t', index=False)

    return 0


# 示例调用
if __name__ == '__main__':
    import pandas as pd

    # 读取示例 TSV 文件
    input = pd.read_csv('phos.tsv', sep='\t')

    # 调用函数
    ptm_sum(input, 'output.tsv')
