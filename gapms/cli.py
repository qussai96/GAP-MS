"""Lightweight command routing; mapping does not import classifier dependencies."""
import sys


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'map':
        from .map_cli import main as map_main
        return map_main(sys.argv[2:])
    from .main import main as workflow_main
    return workflow_main()


if __name__ == '__main__':
    raise SystemExit(main())
