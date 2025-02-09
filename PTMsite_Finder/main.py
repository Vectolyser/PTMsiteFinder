import tkinter as tk
import logging
import sys

from ptm_sitefinder.process_sequence import converting

from tkinter import filedialog, messagebox, ttk
from ptm_sitefinder import UserStoppedException
import threading

# 配置日志
logging.basicConfig(filename='error_log.txt',  # 设置日志文件名
                    filemode='a',  # 追加模式
                    level=logging.ERROR,  # 设置日志级别
                    format='%(asctime)s - %(levelname)s - %(message)s')  # 设置日志格式


def browse_file(entry):
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)


def browse_save_file(entry):
    filename = filedialog.asksaveasfilename(defaultextension=".tsv",
                                            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")])
    entry.delete(0, tk.END)
    entry.insert(0, filename)


def main():

    def start_conversion():
        num = mod_var.get()
        input_tsv_file = tsv_entry.get()
        input_quant_file = quant_entry.get()
        input_fasta_file = fasta_entry.get()
        output_name = output_entry.get()
        filter_threshold = filter_entry.get()

        if not all([num, input_tsv_file, input_quant_file, input_fasta_file, output_name]):
            messagebox.showerror("Error", "All fields are required!")
            return

        # 创建新线程来运行转换函数
        conversion_thread = threading.Thread(target=run_conversion, args=(num, input_tsv_file,
                                                                          input_quant_file, input_fasta_file,
                                                                          output_name, filter_threshold))
        conversion_thread.start()

    def run_conversion(num, input_tsv_file, input_quant_file, input_fasta_file, output_name, filter_threshold):
        try:
            duration = converting(num, input_tsv_file, input_quant_file, input_fasta_file, output_name,
                                  ptm_filter=float(filter_threshold), root=root, progress_bar=progress_bar,
                                  status_label=status_label, start_button=start_button)
            messagebox.showinfo("Success", f"Processing completed!\nDuration: {duration:.4f} seconds")
        except UserStoppedException:
            pass
        except Exception as e:
            logging.exception("An error occurred during the conversion process")
            messagebox.showerror("Error", f"An error occurred: {e}\nPlease check the log for more details.")

        finally:
            start_button.config(text="Start", command=start_conversion)
            if messagebox.askyesno("Exit", "Do you want to close the program?"):
                root.quit()

    def stop_program():
        sys.exit()

    root = tk.Tk()
    root.title("PTMsite Finder")

    tk.Label(root, text="Select Modification:").grid(row=0, column=0, padx=10, pady=5)

    mod_var = tk.StringVar()

    # 创建OptionMenu
    options = ['STY(UniMod:21)', 'K(UniMod:121)', 'K(UniMod:1)', 'N(UniMod:7)', 'K(UniMod:92)']
    mod_var.set(options[0])  # 设置默认值
    mod_menu = tk.OptionMenu(root, mod_var, *options)
    mod_menu.grid(row=1, column=1, padx=10, pady=5, sticky='w')

    # 创建Entry
    mod_entry = tk.Entry(root, textvariable=mod_var, width=15)
    mod_entry.grid(row=0, column=1, padx=10, pady=5, sticky='w')

    mod_entry = tk.Entry(root, textvariable=mod_var, width=15, )
    mod_entry.grid(row=0, column=1, padx=10, pady=5, sticky='w')

    tk.Label(root, text="pr_matrix file:").grid(row=2, column=0, padx=10, pady=5)
    tsv_entry = tk.Entry(root, width=30)
    tsv_entry.grid(row=2, column=1, padx=10, pady=5, sticky='w')
    tk.Button(root, text="Browse", command=lambda: browse_file(tsv_entry)).grid(row=2, column=2, padx=10, pady=5, )

    tk.Label(root, text="main_report file:").grid(row=3, column=0, padx=10, pady=5)
    quant_entry = tk.Entry(root, width=30)
    quant_entry.grid(row=3, column=1, padx=10, pady=5, sticky='w')
    tk.Button(root, text="Browse", command=lambda: browse_file(quant_entry)).grid(row=3, column=2, padx=10, pady=5)

    tk.Label(root, text="Fasta file:").grid(row=4, column=0, padx=10, pady=5)
    fasta_entry = tk.Entry(root, width=30)
    fasta_entry.grid(row=4, column=1, padx=10, pady=5, sticky='w')
    tk.Button(root, text="Browse", command=lambda: browse_file(fasta_entry)).grid(row=4, column=2, padx=10, pady=5)

    tk.Label(root, text="Output file path:").grid(row=5, column=0, padx=10, pady=5)

    tk.Button(root, text="Browse", command=lambda: browse_save_file(output_entry)).grid(row=5, column=2, padx=10,
                                                                                        pady=5)
    output_entry = tk.Entry(root, width=30)
    output_entry.grid(row=5, column=1, padx=10, pady=5, sticky='w')

    tk.Label(root, text="filter threshold:").grid(row=6, column=0, pady=10)
    filter_entry = tk.Entry(root, width=5)
    filter_entry.grid(row=6, column=1, padx=10, pady=5, sticky='w')
    filter_entry.insert(0, "0.75")

    start_button = tk.Button(root, text="Start", command=start_conversion)
    start_button.grid(row=8, column=0, columnspan=3, pady=20)

    progress_bar = ttk.Progressbar(root, orient='horizontal', length=250, mode='determinate')
    progress_bar.grid(row=7, column=1, pady=10)

    status_label = tk.Label(root, text="")
    status_label.grid(row=7, column=2, padx=10, pady=5)

    # 做一个停止运行程序的按钮
    # stop_button = tk.Button(root, text="Stop", command=stop_program           )
    # stop_button.grid(row=8, column=2, padx=10, pady=5)

    root.mainloop()


if __name__ == '__main__':
    main()
