#scripts/yes_no.py

def yes_no(prompt):
    """
    Prompt the user with a yes/no question and return True for 'yes' and False for 'no'.
    
    Args:
        prompt (str): The question to present to the user.
    Returns:
        bool: True if the user answers 'yes', False if 'no'.
    """
    while True:
        choice = input(f"{prompt} (y/n): ").strip().lower()
        if choice in ('y', 'yes'):
            return True
        elif choice in ('n', 'no'):
            return False
        else:
            print("Invalid input. Please enter 'y' for yes or 'n' for no.")