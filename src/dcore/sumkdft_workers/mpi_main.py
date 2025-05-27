import argparse
import sys
import importlib

if __name__ == '__main__':

    try:
        parser = argparse.ArgumentParser(
            description='Internal program for launching SumkDFT.')
        parser.add_argument('runner_cls')
        parser.add_argument('model_hdf5_file')
        parser.add_argument('input_file')
        parser.add_argument('output_file')
        args = parser.parse_args()

        # Runner class
        m_cls = args.runner_cls.split('.')
        if len(m_cls) == 1:
            m = importlib.import_module('dcore.sumkdft_workers.all_workers')
            cls = getattr(m, args.runner_cls)
        else:
            m = importlib.import_module('.'.join(m_cls[:-1]))
            cls = getattr(m, m_cls[-1])
        runner = cls(args.model_hdf5_file, args.input_file, args.output_file)
        runner.run()
    except ImportError as e:
        print("Import error:", e)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    except Exception as e:
        print("Unexpected error:", e)
        import traceback
        traceback.print_exc()
        sys.exit(1)