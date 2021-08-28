import pandas as pd
from post_extraction_analysis_helpers import validation_function


def post_analyse(input_file, output_file):
    df = pd.read_csv(input_file)
    try:
        df['eval'] = df.apply(validation_function, axis=1)
    except ValueError:
        df['eval'] = True
    df.to_csv(output_file,
              sep=',',
              index=False)
