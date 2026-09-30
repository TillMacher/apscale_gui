import subprocess
import sys
import os

def main():
    # Determine theme
    theme = "light"  # default
    if len(sys.argv) > 1:
        if sys.argv[1].lower() in ["dark", "light"]:
            theme = sys.argv[1].lower()

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_dir, 'apscale_gui.py')

        print(f"Starting APSCALE-GUI in {theme} mode. Press CTRL + C to exit.")
        subprocess.run([
            'streamlit', 'run', script_path,
            '--theme.base', theme,
            '--server.address', 'localhost'
        ])

    except KeyboardInterrupt:
        print("Exiting...")
        sys.exit()

if __name__ == "__main__":
    main()