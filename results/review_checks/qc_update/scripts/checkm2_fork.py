"""Wrapper: roda o CheckM2 com multiprocessing 'fork' (padrão no Linux).
No macOS o padrão 'spawn' falha ao serializar Predictor.__set_up_prodigal_thread (AttributeError)."""
import multiprocessing, runpy, sys, shutil
multiprocessing.set_start_method('fork', force=True)
sys.argv[0] = shutil.which('checkm2')
runpy.run_path(sys.argv[0], run_name='__main__')
