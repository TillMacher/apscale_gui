import subprocess, sys, os

def main():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_dir, 'apscale_gui.py')

        print('Press CTRL + C to close APSCALE-GUI!')
        subprocess.run(['streamlit', 'run', script_path, '--theme.base', 'dark', '--server.address', 'localhost'])

    except KeyboardInterrupt:
        print("Exiting...")
        sys.exit()

if __name__ == "__main__":
    main()