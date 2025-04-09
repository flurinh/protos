import torch
import os
import numpy as np
import h5py


def get_embedding(model, tokenizer, device, seq_dict,
                  get_attention=True, rep_dir='models/', save_attention=True,
                  min_seq_len=10, max_seq_len=2000,
                  shift_left=0, shift_right = -1):

    attention_dir = rep_dir + 'attention/'

    embeddings = []

    for seq_idx, (seq_id, seq) in enumerate(seq_dict.items()):
        seq_len = len(seq)
        model_input = [list(seq)]
        if (len(seq) < min_seq_len) or (len(seq) > max_seq_len):
            print("sequence too long!")
        elif len(seq) == 0:
            print("sequence length is ZERO")
        else:
            # add_special_tokens adds extra token at the end of each sequence
            ids = tokenizer.batch_encode_plus(model_input,
                                              add_special_tokens=True,
                                              padding=True,
                                              is_split_into_words=True,
                                              return_tensors="pt")['input_ids'].to(device)
            try:
                with torch.no_grad():
                    # returns: ( batch-size x max_seq_len_in_minibatch x embedding_dim )

                    attn_filename = attention_dir + seq_id + '.h5'
                    if os.path.isfile(attn_filename):
                        get_attns = False
                        save_attns = False
                    else:
                        get_attns = get_attention
                        save_attns = save_attention

                    output = model(input_ids=ids, output_attentions=get_attns)
                    embedding = output.last_hidden_state.detach().cpu().numpy()[shift_left:shift_right]
                    embeddings.append([seq_id, embedding])
                    if save_attns:
                        print("saving attention weights to", attn_filename)
                        attns = [output.attentions[i].detach().cpu().numpy() for i in range(len(output.attentions))]
                        attns = np.stack(attns, axis=1)[0]

                        # open a file, where you ant to store the data
                        h5f = h5py.File(attn_filename, 'w')
                        h5f.create_dataset('attention', data=attns)
                        h5f.create_dataset('embedding', data=embedding)
                        h5f.close()

            except RuntimeError:
                print("RuntimeError during embedding for {} (L={})".format(seq_id, seq_len))
                continue
    return embeddings