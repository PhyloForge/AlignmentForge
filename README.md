# AlignmentForge — phylogenomic alignment & trimming studio

AlignmentForge is an interactive visual software tool for exploring, filtering, trimming, and curating multilocus phylogenomic alignments.

You can inspect sequences, filter them by quality, trim ends, and easily select the best data for your phylogenetic analysis directly in the browser or on your desktop.

**[▶ Try it in your browser](https://phyloforge.github.io/AlignmentForge/?run=example_data)** —
no installation, with example data.

---

## Running it

There are two ways to run AlignmentForge. They are the same application. Use the one you prefer.

| | Use when |
|---|---|
| [1. In a browser](#1-in-a-browser) | Simplest. Nothing to install, works on any OS |
| [2. As a desktop app](#2-as-a-desktop-app) | You want full local filesystem access, offline capability, and maximum performance |

### 1. In a browser

Open **<https://phyloforge.github.io/AlignmentForge/>**. Click the folder icon to open a directory, then select a folder on your computer that contains alignment files.

AlignmentForge uploads nothing. The page reads the folder locally on your machine through the browser's file picker. Your data stays entirely on your device.

### 2. As a desktop app

Download the installer for your operating system (macOS, Windows, or Linux) from the [latest release](https://github.com/PhyloForge/AlignmentForge/releases/latest).

The desktop app gives you complete, unrestricted access to your local filesystem. It can seamlessly read massive folders and write your filtered output datasets back to disk without prompting you for browser permissions.

---

## Try it with the example datasets

You do not need your own alignments to see how AlignmentForge works. This repository ships with several real datasets so you can try out the software right away.

### Exon-only alignments

`example_data/exons/` — Ten exon-only phylip alignment files for frog phylogenomics.

**[▶ Open it live](https://phyloforge.github.io/AlignmentForge/?run=example_data/exons)**

| Get it | How |
|---|---|
| [Browse it on GitHub](https://github.com/PhyloForge/AlignmentForge/tree/main/example_data/exons) | See the sample files |
| [Download the whole repository](https://github.com/PhyloForge/AlignmentForge/archive/refs/heads/main.zip) | `example_data/exons/` is inside it |

### UCE alignments

`example_data/uces/` — Ten Ultraconserved Element (UCE) alignments for testing.

**[▶ Open it live](https://phyloforge.github.io/AlignmentForge/?run=example_data/uces)**

| Get it | How |
|---|---|
| [Browse it on GitHub](https://github.com/PhyloForge/AlignmentForge/tree/main/example_data/uces) | See the sample files |
| [Download the whole repository](https://github.com/PhyloForge/AlignmentForge/archive/refs/heads/main.zip) | `example_data/uces/` is inside it |

### All markers combined

`example_data/all_markers/` — A mix of different marker types in a single folder. 

**[▶ Open it live](https://phyloforge.github.io/AlignmentForge/?run=example_data/all_markers)**

| Get it | How |
|---|---|
| [Browse it on GitHub](https://github.com/PhyloForge/AlignmentForge/tree/main/example_data/all_markers) | See the sample files |
| [Download the whole repository](https://github.com/PhyloForge/AlignmentForge/archive/refs/heads/main.zip) | `example_data/all_markers/` is inside it |

### Loading an example on your machine

1. Download or clone this repository:
```bash
git clone https://github.com/PhyloForge/AlignmentForge.git
```
2. Open AlignmentForge (either the web or desktop version).
3. Click to open a folder.
4. Select `AlignmentForge/example_data/exons` or one of the other example folders.

---

## What to load into AlignmentForge

Load any directory on your computer containing alignment files. 

AlignmentForge currently natively supports standard formats. Once a directory is selected, the application will scan it, calculate summary statistics, and present a visual catalog of all your loci.

You can then apply filters, trim alignments, inspect the alignment matrix view, and export your curated dataset!

---

## Development

If you want to build the software from source:

1. Install **Node.js** (v20+) and **Rust**.
2. Clone this repository.
3. Run `npm install` to install dependencies.
4. Run `npm run dev` to start the browser development server.
5. Run `npm run tauri dev` to start the desktop development application.
