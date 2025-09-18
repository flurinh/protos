from huggingface_hub import login
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

# Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
# login()

protein = ESMProtein(sequence="AAAAA")
client = ESMC.from_pretrained("esmc_300m").to("cuda") # or "cpu"
protein_tensor = client.encode(protein)
logits_output = client.logits(
   protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
)
print(logits_output.logits, logits_output.embeddings)
print(logits_output.embeddings.shape)