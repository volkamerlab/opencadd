def _dic_to_df(dic):
    block_codes = []
    frame_code_categories = []
    frame_code_keywords = []
    data_name_categories = []
    data_name_keywords = []
    data_values = []
    value_dimensions = []

    for block_code, data_block_dic in dic.items():
        for frame_code_category, save_frame_category_dic in data_block_dic.items():
            for frame_code_keyword, save_frame_keyword_dic in save_frame_category_dic.items():
                for data_name_category, data_category_dic in save_frame_keyword_dic.items():
                    for data_name_keyword, data_value in data_category_dic.items():
                        block_codes.append(block_code)
                        frame_code_categories.append(frame_code_category)
                        frame_code_keywords.append(frame_code_keyword)
                        data_name_categories.append(data_name_category)
                        data_name_keywords.append(data_name_keyword)
                        data_values.append(data_value)
                        value_dimensions.append(value_dimension)

    data = dict(
        block_code=block_codes,
        frame_code_category=frame_code_categories,
        frame_code_keyword=frame_code_keywords,
        data_name_category=data_name_categories,
        data_name_keyword=data_name_keywords,
        data_value=data_values,
        value_dimension=value_dimensions,
    )
    return pd.DataFrame(data)