from pyparsing import Word, alphas, infixNotation, opAssoc, oneOf, nums


class LogicParser:
    def __init__(self):
        self.identifier = Word(nums + alphas + "_" + "/" + "-" + "*" + ".") # + " "
        self.operators = oneOf("and or not")

        self.expr = infixNotation(
            self.identifier,
            [
                ("not", 1, opAssoc.RIGHT),  # ?????
                ("and", 2, opAssoc.LEFT),  # ?????????
                ("or", 2, opAssoc.LEFT),  # ?????????
            ],
        )

    def parse(self, logic_string):
        return self.expr.parseString(logic_string).asList()

    def extract_elements(self, parsed_list):
        flat_list = []
        stack = list(parsed_list)
        while stack:
            item = stack.pop(0)
            if isinstance(item, list):
                stack = item + stack
            else:
                flat_list.append(item)

        seen = set()
        elements = []
        for item in flat_list:
            if item not in {"and", "or", "not"} and item not in seen:
                seen.add(item)
                elements.append(item)
        return elements

    def generate_mapping(self, parsed_list):
        return {elem: f"E{idx}" for idx, elem in enumerate(self.extract_elements(parsed_list), start=1)}

    def extract_structure(self, parsed_list, mapping):
        replaced = []
        for item in parsed_list:
            if isinstance(item, list):  # ?????????
                replaced.append(self.extract_structure(item, mapping))
            elif item in mapping:  # ?????
                replaced.append(mapping[item])
            else:  # ?????
                replaced.append(item)
        return replaced

    # def extract_elements_(self, parsed_list):
    #     elements = set()
    #     for item in parsed_list:
    #         if isinstance(item, list):  # ????????????
    #             elements.update(self.extract_elements_(item))
    #         elif item not in {"and", "or", "not"}:  # ??????????
    #             elements.add(item)
    #     # elements = list(self.generate_mapping(parsed_list).keys())
    #     return elements
    #
    # def generate_mapping_(self, parsed_list):
    #     elements = set()
    #
    #     def extract_elements(sublist):
    #         for item in sublist:
    #             if isinstance(item, list):  # ?????????
    #                 extract_elements(item)
    #             elif item not in {"and", "or", "not"}:  # ????
    #                 elements.add(item)
    #
    #     extract_elements(parsed_list)
    #     return {elem: f"E{idx}" for idx, elem in enumerate(sorted(elements), start=1)}


    def reconstruct(self, parsed_list):
        reconstructed = []
        for item in parsed_list:
            if isinstance(item, list):  # ?????????
                reconstructed.append(f"({self.reconstruct(item)})")
            else:
                reconstructed.append(item)
        return " ".join(reconstructed)


if __name__ == "__main__":
    logic_string = "(EGFR_Exon19del or EGFR_L858R or EGFR_T790M) and (CTNNB1_unspecified or PIK3CA_unspecified or RB1_unspecified or TP53_unspecified)"

    parser = LogicParser()

    parsed_result1 = parser.parse(logic_string)
    print("????:", parsed_result1)

    elements1 = parser.extract_elements(parsed_result1)
    print("????:", elements1)

    mapping1 = parser.generate_mapping(parsed_result1)
    structure1 = parser.extract_structure(parsed_result1, mapping1)
    print("????:", structure1)

    str1 = parser.reconstruct(parsed_result1)
    str1 = str1[1:-1]
    print('ok')

