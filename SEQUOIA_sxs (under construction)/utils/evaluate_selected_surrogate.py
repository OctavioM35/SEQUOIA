from config import surrogate_model
from utils.surrogates import surrogates

def evaluate_selected_surrogate(surrogate_model):
    for precessing_or_not in surrogates.keys():
        for suported_surrogate_particular in surrogates[precessing_or_not]:
            if surrogate_model.lower() in suported_surrogate_particular.lower():
                print('\nSelected surrogate accepted. \nSelected surrogate model: ' , suported_surrogate_particular + "\n")
                if precessing_or_not == 'non-precessing':
                    return False
                else:
                    return True
    message = surrogate_model + " surrogate model not supported by SEQUOIA.\nSupported models:\n"

    for category, models in surrogates.items():
        message += f"  {category}:\n"
        for model in models:
            message += f"    - {model}\n"

    raise ValueError(message)