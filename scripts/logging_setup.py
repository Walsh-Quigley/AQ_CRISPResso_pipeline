# scripts/logging.py
import logging
import os
from datetime import datetime

def setup_logging():
    """Configure logging to write to both file and console."""
    # Create logs directory if it doesn't exist
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Generate log filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"CRISPResso_run_{timestamp}.log")
    
    # Configure logging format
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # This still prints to console
        ]
    )
    
    logging.info(f"Log file created: {log_file}")
    return log_file