#!/usr/bin/env python3
import subprocess
import os
import sys

def run_matlab_script():
    """Execute the RNA-seq MATLAB algorithm"""
    
    # Change to the RNA-seq directory
    rna_seq_dir = "/home/marco/projects/TCC_Inference_Methods/rna-seq-main"
    os.chdir(rna_seq_dir)
    
    # MATLAB command to run the main algorithm
    matlab_cmd = [
        "/home/marco/matlab/bin/matlab", 
        "-nodisplay", 
        "-nosplash", 
        "-nodesktop",
        "-r", 
        "cd('code'); run('GCSTI_main code.m'); exit;"
    ]
    
    try:
        print("Executando algoritmo RNA-seq...")
        print("Diretório:", os.getcwd())
        
        # Run MATLAB
        result = subprocess.run(matlab_cmd, 
                              capture_output=True, 
                              text=True, 
                              timeout=300)
        
        if result.returncode == 0:
            print("✓ Algoritmo executado com sucesso!")
            print("\nSaída:")
            print(result.stdout)
        else:
            print("✗ Erro na execução:")
            print(result.stderr)
            
    except subprocess.TimeoutExpired:
        print("✗ Timeout - execução demorou mais que 5 minutos")
    except FileNotFoundError:
        print("✗ MATLAB não encontrado. Certifique-se que está instalado e no PATH")
    except Exception as e:
        print(f"✗ Erro: {e}")

if __name__ == "__main__":
    run_matlab_script()