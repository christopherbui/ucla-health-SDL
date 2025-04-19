import sys
import csv
from pathlib import Path
from typing import Union


def txt_to_csv(source_file: Union[str, Path], output_file: Union[str, Path], delim: str ='\t'):
    """Convert a delimited .txt file to .csv
    
    Parameters:
        source_file: text file
        output_file: csv file
        delim: delimeter
    
    Return:
        None
    """
    with open(source_file, 'r', encoding='utf-8') as txt_file:
        reader = csv.reader(txt_file, delimiter=delim)
        with open(output_file, 'w', newline='', encoding='utf-8') as csv_file:
            writer = csv.writer(csv_file)
            for row in reader:
                writer.writerow(row)
    print(f"[✓✓✓] Converted '{source_file}' to '{output_file}' as CSV.\n")


def main():
    # make sure enough arguments are passed
    if len(sys.argv) < 4:
        print("Usage: python converter.py <conversion_function> <source_file> <output_file>\n")
        return
    
    conversion_type = sys.argv[1]
    source_file = sys.argv[2]
    output_file = sys.argv[3]

    print("Attempting Conversion...\n")

    match conversion_type:
        case 'txt_to_csv':
            txt_to_csv(source_file, output_file)
        case _:
            print(f"[!!!] Converstion type '{conversion_type}' not supported!\n")


if __name__ == '__main__':
    main()