import pandas as pd
import os
from tkinter import messagebox
from Bio import SeqIO
from PTMsite_Finder import converting_rules, timing_decorator, UserStoppedException
from PTMsite_Finder.process_PTM import ptm_sum
import concurrent.futures

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('max_colwidth', 200)

running_flag = True
stop_confirmed = False


def stop_conversion():
    global running_flag
    running_flag = False


@timing_decorator
def converting(num, input_tsv_file, input_quant_file, input_fasta_file, output_name, ptm_filter=0.75, root=None,
               progress_bar=None, status_label=None, start_button=None):
    progress_bar = progress_bar

    start_button.config(text="Stop", command=stop_conversion)

    def update_progress(value, message=""):
        progress_bar['value'] = value
        if status_label is not None:
            status_label.config(text=message)
        root.update_idletasks()

    def check_state():
        global running_flag
        if not running_flag:
            running_flag = True
            messagebox.showinfo("Info", "Conversion stopped by user")
            update_progress(0, "stopped")
            raise UserStoppedException

    update_progress(0, "start processing...")

    pr_matrix = pd.read_csv(input_tsv_file, sep='\t')

    check_state()

    tsv = pd.read_csv(input_quant_file, sep='\t')

    # 保留Modified.Amino.Acid列有值的且PTM.Site.Confidence > 0.75的行
    tsv = tsv[tsv['PTM.Site.Confidence'] >= ptm_filter]

    update_progress(20, "loading fasta...")

    fasta_dict = {}
    for protein in SeqIO.parse(open(input_fasta_file), 'fasta'):

        check_state()

        protein_id_parts = protein.id.split('|')
        if len(protein_id_parts) > 1:  # 检查分割后的列表长度
            fasta_dict[protein_id_parts[1]] = str(protein.seq)
        else:
            print(f"Warning: Protein ID '{protein.id}' does not contain '|' and will be skipped.")

    update_progress(30, "merging data...")

    tsv = tsv[['Protein.Group', 'Precursor.Id', 'PTM.Site.Confidence']].copy()
    ori_df = pd.merge(pr_matrix, tsv, on=['Protein.Group', 'Precursor.Id'], how='inner')

    check_state()

    df = ori_df[['Protein.Group', 'Stripped.Sequence', 'Modified.Sequence']].copy()

    update_progress(40, "merging data...")

    df.loc[:, 'Modified.Sequence'] = df['Modified.Sequence'].apply(lambda x: converting_rules(x, num=num))
    df_dict = df.to_dict('records')

    update_progress(50, "analysing data...")

    # 新增：将数据处理分成函数并实行多线程
    def process_data(data):
        Id = data['Protein.Group'].split(';')[0]
        seq = data['Stripped.Sequence']
        Mod_seq = data['Modified.Sequence']

        if Id in fasta_dict:

            check_state()

            pro_seq = fasta_dict[Id]
            start = pro_seq.find(seq)
            if start == -1:
                if 'I' in seq:
                    seq = seq.replace('I', 'L')
                if 'L' in seq:
                    seq = seq.replace('L', 'I')

                start = pro_seq.find(seq)

            data['Modified.Amino.Acid'] = [seq[i] for i, char in enumerate(Mod_seq) if char.isdigit()]
            data['Peptide.Position'] = [int(i) + 1 for i, char in enumerate(Mod_seq) if char.isdigit()]

            data['Protein.Position'] = [int(i) + int(start) for i in data['Peptide.Position']]
            sequence_windows = [pro_seq[max(0, i - 8):i] + '*' + pro_seq[i: i + 7] for i in data['Protein.Position']]
            data['Sequence.Window'] = sequence_windows
            data['Start'] = start + 1
            data['End'] = int(start + 1) + len(seq)

            return data  # 返回数据以便后续处理
        else:
            data.update({
                'Modified.Amino.Acid': None,
                'Peptide.Position': None,
                'Protein.Position': None,
                'Sequence.Window': None,
                'Start': None,
                'End': None
            })
            return data

    with concurrent.futures.ThreadPoolExecutor() as executor:  # 新增：创建线程池
        df_dict = list(executor.map(process_data, df.to_dict('records')))  # 将处理数据的函数应用到每条记录上

    count_list = []
    for data in df_dict:

        check_state()

        if data['Modified.Amino.Acid'] is not None:
            count_list.extend(
                [data['Protein.Group'].split(';')[0] + '_' + str(pos) for pos in data['Protein.Position']])

    update_progress(70, "preparing for output...")

    out_df = pd.DataFrame(df_dict)
    out_df['Start'] = pd.to_numeric(out_df['Start'], errors='coerce').astype('Int64')
    out_df['End'] = pd.to_numeric(out_df['End'], errors='coerce').astype('Int64')
    out_df['Sequence.Window'] = out_df['Sequence.Window'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)
    out_df['Modified.Amino.Acid'] = out_df['Modified.Amino.Acid'].apply(
        lambda x: ','.join(x) if isinstance(x, list) else x)
    out_df['Peptide.Position'] = out_df['Peptide.Position'].apply(
        lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)
    out_df['Protein.Position'] = out_df['Protein.Position'].apply(
        lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)

    out_df.drop(columns=['Modified.Sequence', 'Protein.Group', 'Stripped.Sequence'], inplace=True)

    check_state()

    concat_df = pd.concat([ori_df, out_df], axis=1)

    check_state()

    concat_df['Protein.Position.Split'] = concat_df['Protein.Position'].apply(lambda x: str(x).split(', '))
    concat_df['Modified.Amino.Acid.Split'] = concat_df['Modified.Amino.Acid'].apply(lambda x: str(x).split(','))
    concat_df['Combined.Split'] = concat_df.apply(
        lambda row: [f"{a}{b}" for a, b in zip(row['Modified.Amino.Acid.Split'], row['Protein.Position.Split'])],
        axis=1)

    concat_df['Protein.Group_Position'] = concat_df.apply(
        lambda row: [f"{row['Protein.Group'].split(';')[0]}_{pos}" for pos in row['Combined.Split']] if row[
            'Modified.Amino.Acid'] else None, axis=1)

    check_state()

    concat_df.sort_values('PTM.Site.Confidence', ascending=False, inplace=True)
    concat_df.reset_index(drop=True, inplace=True)
    concat_df.drop_duplicates(subset=['Precursor.Id', 'Protein.Group'], inplace=True, keep='first')

    # 更新状态信息
    update_progress(80, "calculating PTM...")

    check_state()

    concat_df['Unique PTM site'] = None
    concat_df.loc[0, 'Unique PTM site'] = len(set(count_list))

    concat_df['Protein.Group_Position'] = concat_df['Protein.Group_Position'].apply(lambda x: ', '.join(x) if
    x else '')
    df_grouped = concat_df.groupby('Protein.Group_Position')['PTM.Site.Confidence'].max().reset_index()

    num_rows = df_grouped[df_grouped['PTM.Site.Confidence'] > 0.75].shape[0]

    concat_df['PTM > 0.75 count'] = None
    concat_df.loc[0, 'PTM > 0.75 count'] = num_rows

    check_state()

    update_progress(90, "saving output...")

    concat_df.drop(columns=['Protein.Position.Split', 'Modified.Amino.Acid.Split', 'Combined.Split'], inplace=True)
    base_name = os.path.splitext(output_name)[0]
    concat_df.to_csv(r'{}'.format(output_name), sep='\t', index=False)
    ptm_sum(concat_df, output_file=f'{base_name}_ptm.tsv')

    update_progress(100, "completed")

    return 0
