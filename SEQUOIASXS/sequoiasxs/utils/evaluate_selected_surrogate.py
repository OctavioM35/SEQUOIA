from utils.surrogates import surrogates

def evaluate_selected_surrogate(surrogate_model):
    """
    Check that the selected surrogate model is supported and determine
    whether it includes spin precession.

    The function searches the registry of supported surrogates
    (utils.surrogates.surrogates), which groups the models into
    'non-precessing' and precessing categories. The comparison is
    case-insensitive and succeeds if the given name is contained in the name
    of a registered model. When a match is found, the selected model is
    printed. If no match is found, a ValueError listing all the supported
    models by category is raised.

    INPUTS:
    -------

    -surrogate_model: Name of the surrogate model selected by the user
    (e.g. "NRSur7dq4"). A partial name is accepted as long as it is
    contained in the name of a supported model.

    OUTPUTS:
    --------

    -False if the matched surrogate belongs to the 'non-precessing' category.

    -True if it belongs to any other category (precessing).

    RAISES:
    -------

    -ValueError: if the surrogate model is not supported by SEQUOIA. The
    message includes the list of supported models.
    """
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